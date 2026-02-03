use std::collections::HashMap;
use std::io::BufRead;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::error::{Result, SdfError};
use crate::molecule::{Molecule, SdfFormat};

/// Detected file format for molecular structure files.
///
/// This enum represents all supported file formats for automatic detection
/// and parsing. It is used by [`detect_format`] and the `parse_auto_*` functions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    /// SDF V2000 format (traditional SDF, up to 999 atoms/bonds)
    SdfV2000,
    /// SDF V3000 format (extended format, no atom/bond limits)
    SdfV3000,
    /// TRIPOS MOL2 format
    Mol2,
    /// XYZ format (coordinates only, no bonds)
    Xyz,
}

impl std::fmt::Display for FileFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FileFormat::SdfV2000 => write!(f, "sdf_v2000"),
            FileFormat::SdfV3000 => write!(f, "sdf_v3000"),
            FileFormat::Mol2 => write!(f, "mol2"),
            FileFormat::Xyz => write!(f, "xyz"),
        }
    }
}

/// SDF V2000 format parser.
pub struct SdfParser<R> {
    reader: R,
    line_number: usize,
}

impl<R: BufRead> SdfParser<R> {
    /// Creates a new parser from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
        }
    }

    /// Reads the next line from the input.
    fn read_line(&mut self, buf: &mut String) -> Result<bool> {
        buf.clear();
        let bytes_read = self.reader.read_line(buf)?;
        if bytes_read > 0 {
            self.line_number += 1;
            // Remove trailing newline
            if buf.ends_with('\n') {
                buf.pop();
                if buf.ends_with('\r') {
                    buf.pop();
                }
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Parses a single molecule from the input.
    /// Returns None if end of file is reached.
    pub fn parse_molecule(&mut self) -> Result<Option<Molecule>> {
        let mut line = String::new();

        // Line 1: Molecule name
        if !self.read_line(&mut line)? {
            return Ok(None);
        }
        let name = line.trim().to_string();

        // Line 2: Program/timestamp line
        if !self.read_line(&mut line)? {
            return Err(SdfError::MissingSection("header".to_string()));
        }
        let program_line = if line.trim().is_empty() {
            None
        } else {
            Some(line.clone())
        };

        // Line 3: Comment
        if !self.read_line(&mut line)? {
            return Err(SdfError::MissingSection("header".to_string()));
        }
        let comment = if line.trim().is_empty() {
            None
        } else {
            Some(line.clone())
        };

        // Line 4: Counts line
        if !self.read_line(&mut line)? {
            return Err(SdfError::MissingSection("counts line".to_string()));
        }
        let (atom_count, bond_count) = self.parse_counts_line(&line)?;

        // Parse atoms
        let mut atoms = Vec::with_capacity(atom_count);
        for i in 0..atom_count {
            if !self.read_line(&mut line)? {
                return Err(SdfError::AtomCountMismatch {
                    expected: atom_count,
                    found: i,
                });
            }
            let atom = self.parse_atom_line(&line, i)?;
            atoms.push(atom);
        }

        // Parse bonds
        let mut bonds = Vec::with_capacity(bond_count);
        for i in 0..bond_count {
            if !self.read_line(&mut line)? {
                return Err(SdfError::BondCountMismatch {
                    expected: bond_count,
                    found: i,
                });
            }
            let bond = self.parse_bond_line(&line, atom_count)?;
            bonds.push(bond);
        }

        // Parse property block until M  END
        let mut properties = HashMap::new();
        loop {
            if !self.read_line(&mut line)? {
                break;
            }
            if line.starts_with("M  END") {
                break;
            }
            // Handle M  CHG (charge) lines
            if line.starts_with("M  CHG") {
                self.parse_charge_line(&line, &mut atoms)?;
            }
            // Handle M  ISO (isotope) lines
            if line.starts_with("M  ISO") {
                self.parse_isotope_line(&line, &mut atoms)?;
            }
        }

        // Parse data block until $$$$ or EOF
        let mut current_property_name: Option<String> = None;
        let mut current_property_value = String::new();

        loop {
            if !self.read_line(&mut line)? {
                break;
            }
            if line.starts_with("$$$$") {
                // Save any pending property
                if let Some(prop_name) = current_property_name.take() {
                    properties.insert(prop_name, current_property_value.trim().to_string());
                }
                break;
            }
            if line.starts_with("> ") || line.starts_with(">  ") {
                // Save previous property if exists
                if let Some(prop_name) = current_property_name.take() {
                    properties.insert(prop_name, current_property_value.trim().to_string());
                }
                current_property_value.clear();

                // Extract property name from angle brackets
                if let Some(start) = line.find('<') {
                    // Find closing '>' after the opening '<'
                    if let Some(end) = line[start + 1..].find('>') {
                        let prop_name = line[start + 1..start + 1 + end].to_string();
                        current_property_name = Some(prop_name);
                    }
                }
            } else if current_property_name.is_some() && !line.is_empty() {
                // Accumulate property value
                if !current_property_value.is_empty() {
                    current_property_value.push('\n');
                }
                current_property_value.push_str(&line);
            }
        }

        Ok(Some(Molecule {
            name,
            program_line,
            comment,
            atoms,
            bonds,
            properties,
            format_version: SdfFormat::V2000,
            stereogroups: Vec::new(),
            sgroups: Vec::new(),
            collections: Vec::new(),
        }))
    }

    /// Parses the counts line (line 4 of the molfile).
    fn parse_counts_line(&self, line: &str) -> Result<(usize, usize)> {
        // V2000 format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
        // Positions: 0-2 = atom count, 3-5 = bond count
        if line.len() < 6 {
            return Err(SdfError::InvalidCountsLine(line.to_string()));
        }

        let atom_count: usize = line[0..3]
            .trim()
            .parse()
            .map_err(|_| SdfError::InvalidCountsLine(line.to_string()))?;

        let bond_count: usize = line[3..6]
            .trim()
            .parse()
            .map_err(|_| SdfError::InvalidCountsLine(line.to_string()))?;

        Ok((atom_count, bond_count))
    }

    /// Parses an atom line.
    fn parse_atom_line(&self, line: &str, index: usize) -> Result<Atom> {
        // V2000 atom line format:
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // Positions: 0-9 x, 10-19 y, 20-29 z, 31-33 symbol, 34-35 mass diff, 36-38 charge

        if line.len() < 34 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Atom line too short: {}", line),
            });
        }

        let x: f64 = line[0..10]
            .trim()
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(line[0..10].to_string()))?;

        let y: f64 = line[10..20]
            .trim()
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(line[10..20].to_string()))?;

        let z: f64 = line[20..30]
            .trim()
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(line[20..30].to_string()))?;

        let element = line[31..34].trim().to_string();

        // Parse optional fields
        let mass_difference: i8 = if line.len() >= 36 {
            line[34..36].trim().parse().unwrap_or(0)
        } else {
            0
        };

        // Charge mapping in V2000: 0=uncharged, 1=+3, 2=+2, 3=+1, 4=doublet radical, 5=-1, 6=-2, 7=-3
        let formal_charge: i8 = if line.len() >= 39 {
            let charge_code: u8 = line[36..39].trim().parse().unwrap_or(0);
            match charge_code {
                0 => 0,
                1 => 3,
                2 => 2,
                3 => 1,
                4 => 0, // doublet radical
                5 => -1,
                6 => -2,
                7 => -3,
                _ => 0,
            }
        } else {
            0
        };

        // Stereo parity
        let stereo_parity: Option<u8> = if line.len() >= 42 {
            let parity: u8 = line[39..42].trim().parse().unwrap_or(0);
            if parity > 0 { Some(parity) } else { None }
        } else {
            None
        };

        // Hydrogen count
        let hydrogen_count: Option<u8> = if line.len() >= 45 {
            let hcount: u8 = line[42..45].trim().parse().unwrap_or(0);
            if hcount > 0 { Some(hcount) } else { None }
        } else {
            None
        };

        // Valence
        let valence: Option<u8> = if line.len() >= 51 {
            let val: u8 = line[48..51].trim().parse().unwrap_or(0);
            if val > 0 { Some(val) } else { None }
        } else {
            None
        };

        Ok(Atom {
            index,
            element,
            x,
            y,
            z,
            formal_charge,
            mass_difference,
            stereo_parity,
            hydrogen_count,
            valence,
            v3000_id: None,
            atom_atom_mapping: None,
            rgroup_label: None,
            radical: None,
        })
    }

    /// Parses a bond line.
    fn parse_bond_line(&self, line: &str, atom_count: usize) -> Result<Bond> {
        // V2000 bond line format:
        // 111222tttsssxxxrrrccc
        // Positions: 0-2 = first atom, 3-5 = second atom, 6-8 = bond type, 9-11 = stereo

        if line.len() < 9 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Bond line too short: {}", line),
            });
        }

        let atom1: usize = line[0..3]
            .trim()
            .parse::<usize>()
            .map_err(|_| SdfError::Parse {
                line: self.line_number,
                message: "Invalid atom1 index".to_string(),
            })?;

        let atom2: usize = line[3..6]
            .trim()
            .parse::<usize>()
            .map_err(|_| SdfError::Parse {
                line: self.line_number,
                message: "Invalid atom2 index".to_string(),
            })?;

        // Convert from 1-based to 0-based indices
        let atom1 = atom1.checked_sub(1).ok_or(SdfError::InvalidAtomIndex {
            index: atom1,
            atom_count,
        })?;

        let atom2 = atom2.checked_sub(1).ok_or(SdfError::InvalidAtomIndex {
            index: atom2,
            atom_count,
        })?;

        // Validate atom indices
        if atom1 >= atom_count {
            return Err(SdfError::InvalidAtomIndex {
                index: atom1 + 1,
                atom_count,
            });
        }
        if atom2 >= atom_count {
            return Err(SdfError::InvalidAtomIndex {
                index: atom2 + 1,
                atom_count,
            });
        }

        let bond_type: u8 = line[6..9].trim().parse().map_err(|_| SdfError::Parse {
            line: self.line_number,
            message: "Invalid bond type".to_string(),
        })?;

        let order = BondOrder::from_sdf(bond_type).ok_or(SdfError::InvalidBondOrder(bond_type))?;

        let stereo = if line.len() >= 12 {
            let stereo_code: u8 = line[9..12].trim().parse().unwrap_or(0);
            BondStereo::from_sdf(stereo_code)
        } else {
            BondStereo::None
        };

        let topology = if line.len() >= 18 {
            let topo: u8 = line[15..18].trim().parse().unwrap_or(0);
            if topo > 0 { Some(topo) } else { None }
        } else {
            None
        };

        Ok(Bond {
            atom1,
            atom2,
            order,
            stereo,
            topology,
            v3000_id: None,
            reacting_center: None,
        })
    }

    /// Parses M  CHG charge lines and updates atoms.
    fn parse_charge_line(&self, line: &str, atoms: &mut [Atom]) -> Result<()> {
        // Format: M  CHG  n   aaa vvv   aaa vvv ...
        // n = number of entries, aaa = atom number (1-based), vvv = charge
        if line.len() < 9 {
            return Ok(());
        }

        let count: usize = line[6..9].trim().parse().unwrap_or(0);
        let mut pos = 9;

        for _ in 0..count {
            if pos + 8 > line.len() {
                break;
            }
            let atom_num: usize = line[pos..pos + 4].trim().parse().unwrap_or(0);
            let charge: i8 = line[pos + 4..pos + 8].trim().parse().unwrap_or(0);

            if atom_num > 0 && atom_num <= atoms.len() {
                atoms[atom_num - 1].formal_charge = charge;
            }
            pos += 8;
        }

        Ok(())
    }

    /// Parses M  ISO isotope lines and updates atoms.
    fn parse_isotope_line(&self, line: &str, atoms: &mut [Atom]) -> Result<()> {
        // Format: M  ISO  n   aaa vvv   aaa vvv ...
        if line.len() < 9 {
            return Ok(());
        }

        let count: usize = line[6..9].trim().parse().unwrap_or(0);
        let mut pos = 9;

        for _ in 0..count {
            if pos + 8 > line.len() {
                break;
            }
            let atom_num: usize = line[pos..pos + 4].trim().parse().unwrap_or(0);
            let mass_diff: i8 = line[pos + 4..pos + 8].trim().parse().unwrap_or(0);

            if atom_num > 0 && atom_num <= atoms.len() {
                atoms[atom_num - 1].mass_difference = mass_diff;
            }
            pos += 8;
        }

        Ok(())
    }
}

/// Iterator over molecules in an SDF file.
pub struct SdfIterator<R> {
    parser: SdfParser<R>,
    finished: bool,
}

impl<R: BufRead> SdfIterator<R> {
    /// Creates a new iterator from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            parser: SdfParser::new(reader),
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for SdfIterator<R> {
    type Item = Result<Molecule>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }

        match self.parser.parse_molecule() {
            Ok(Some(mol)) => Some(Ok(mol)),
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

/// Parses a single molecule from an SDF string.
pub fn parse_sdf_string(content: &str) -> Result<Molecule> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let mut parser = SdfParser::new(reader);

    parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
}

/// Parses all molecules from an SDF string.
pub fn parse_sdf_string_multi(content: &str) -> Result<Vec<Molecule>> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let iter = SdfIterator::new(reader);

    iter.collect()
}

/// Parses a single molecule from an SDF file.
pub fn parse_sdf_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut parser = SdfParser::new(reader);

    parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
}

/// Parses all molecules from an SDF file.
pub fn parse_sdf_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let iter = SdfIterator::new(reader);

    iter.collect()
}

/// Returns an iterator over molecules in an SDF file.
pub fn iter_sdf_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<SdfIterator<std::io::BufReader<std::fs::File>>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    Ok(SdfIterator::new(reader))
}

/// Detects the SDF format version from file content.
///
/// Reads the first few lines to check for V3000 indicators.
pub fn detect_sdf_format(content: &str) -> SdfFormat {
    // V3000 files have "V3000" in the counts line (line 4)
    let lines: Vec<&str> = content.lines().take(5).collect();
    if lines.len() >= 4 && lines[3].contains("V3000") {
        SdfFormat::V3000
    } else {
        SdfFormat::V2000
    }
}

/// Parses an SDF string with automatic format detection.
///
/// This function automatically detects whether the content is V2000 or V3000
/// format and uses the appropriate parser.
pub fn parse_sdf_auto_string(content: &str) -> Result<Molecule> {
    match detect_sdf_format(content) {
        SdfFormat::V2000 => parse_sdf_string(content),
        SdfFormat::V3000 => super::sdf_v3000::parse_sdf_v3000_string(content),
    }
}

/// Parses multiple molecules from an SDF string with automatic format detection.
pub fn parse_sdf_auto_string_multi(content: &str) -> Result<Vec<Molecule>> {
    match detect_sdf_format(content) {
        SdfFormat::V2000 => parse_sdf_string_multi(content),
        SdfFormat::V3000 => super::sdf_v3000::parse_sdf_v3000_string_multi(content),
    }
}

/// Parses an SDF file with automatic format detection.
pub fn parse_sdf_auto_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    // Read the first few lines to detect format
    let content = std::fs::read_to_string(&path)?;
    match detect_sdf_format(&content) {
        SdfFormat::V2000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            let mut parser = SdfParser::new(reader);
            parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
        }
        SdfFormat::V3000 => super::sdf_v3000::parse_sdf_v3000_file(path),
    }
}

/// Parses all molecules from an SDF file with automatic format detection.
pub fn parse_sdf_auto_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    let content = std::fs::read_to_string(&path)?;
    match detect_sdf_format(&content) {
        SdfFormat::V2000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            let iter = SdfIterator::new(reader);
            iter.collect()
        }
        SdfFormat::V3000 => super::sdf_v3000::parse_sdf_v3000_file_multi(path),
    }
}

// ============================================================================
// Unified Auto-Detection (SDF V2000, V3000, MOL2)
// ============================================================================

/// Detects the format of a molecular structure file from its content.
///
/// This function examines the content to determine whether it is:
/// - MOL2 format (contains `@<TRIPOS>` marker)
/// - XYZ format (first line is integer, third line has element + 3 floats)
/// - SDF V3000 format (contains `V3000` on the counts line)
/// - SDF V2000 format (default)
///
/// Detection order: MOL2 → XYZ → V3000 → V2000 (default)
///
/// # Example
///
/// ```rust
/// use sdfrust::detect_format;
/// use sdfrust::FileFormat;
///
/// let mol2_content = "@<TRIPOS>MOLECULE\ntest\n";
/// assert_eq!(detect_format(mol2_content), FileFormat::Mol2);
///
/// let xyz_content = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n";
/// assert_eq!(detect_format(xyz_content), FileFormat::Xyz);
///
/// let v2000_content = "test\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n";
/// assert_eq!(detect_format(v2000_content), FileFormat::SdfV2000);
/// ```
pub fn detect_format(content: &str) -> FileFormat {
    // Check first ~100 lines for @<TRIPOS> marker (MOL2)
    for line in content.lines().take(100) {
        if line.starts_with("@<TRIPOS>") {
            return FileFormat::Mol2;
        }
    }

    // Check for XYZ format:
    // - First line should be an integer (atom count)
    // - Third line should be: element x y z (element + 3 floats)
    let lines: Vec<&str> = content.lines().take(5).collect();
    if lines.len() >= 3 {
        let first_line_is_int = lines[0].trim().parse::<usize>().is_ok();
        if first_line_is_int {
            // Check if line 3 looks like an XYZ atom line
            let parts: Vec<&str> = lines[2].split_whitespace().collect();
            if parts.len() >= 4 {
                // First part should be element (letters or small number)
                // Next three should be parseable as floats
                let looks_like_element =
                    parts[0].chars().all(|c| c.is_alphabetic()) || parts[0].parse::<u8>().is_ok();
                let has_three_coords = parts[1].parse::<f64>().is_ok()
                    && parts[2].parse::<f64>().is_ok()
                    && parts[3].parse::<f64>().is_ok();

                if looks_like_element && has_three_coords {
                    return FileFormat::Xyz;
                }
            }
        }
    }

    // Check line 4 for V3000 marker
    if lines.len() >= 4 && lines[3].contains("V3000") {
        return FileFormat::SdfV3000;
    }

    // Default to V2000
    FileFormat::SdfV2000
}

/// Parses a single molecule from a string with automatic format detection.
///
/// This function automatically detects whether the content is SDF V2000, V3000,
/// MOL2, or XYZ format and uses the appropriate parser.
///
/// # Example
///
/// ```rust,ignore
/// use sdfrust::parse_auto_string;
///
/// let mol = parse_auto_string(sdf_content).unwrap();
/// assert_eq!(mol.name, "methane");
/// ```
pub fn parse_auto_string(content: &str) -> Result<Molecule> {
    match detect_format(content) {
        FileFormat::SdfV2000 => parse_sdf_string(content),
        FileFormat::SdfV3000 => super::sdf_v3000::parse_sdf_v3000_string(content),
        FileFormat::Mol2 => super::mol2::parse_mol2_string(content),
        FileFormat::Xyz => super::xyz::parse_xyz_string(content),
    }
}

/// Parses multiple molecules from a string with automatic format detection.
///
/// This function automatically detects whether the content is SDF V2000, V3000,
/// MOL2, or XYZ format and uses the appropriate parser.
pub fn parse_auto_string_multi(content: &str) -> Result<Vec<Molecule>> {
    match detect_format(content) {
        FileFormat::SdfV2000 => parse_sdf_string_multi(content),
        FileFormat::SdfV3000 => super::sdf_v3000::parse_sdf_v3000_string_multi(content),
        FileFormat::Mol2 => super::mol2::parse_mol2_string_multi(content),
        FileFormat::Xyz => super::xyz::parse_xyz_string_multi(content),
    }
}

/// Parses a single molecule from a file with automatic format detection.
///
/// This function reads the file, detects whether it is SDF V2000, V3000,
/// MOL2, or XYZ format, and uses the appropriate parser.
///
/// # Example
///
/// ```rust,ignore
/// use sdfrust::parse_auto_file;
///
/// // Works with any supported format
/// let mol = parse_auto_file("molecule.sdf")?;  // SDF V2000 or V3000
/// let mol = parse_auto_file("molecule.mol2")?; // MOL2
/// let mol = parse_auto_file("molecule.xyz")?;  // XYZ
/// ```
pub fn parse_auto_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    let content = std::fs::read_to_string(&path)?;
    match detect_format(&content) {
        FileFormat::SdfV2000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            let mut parser = SdfParser::new(reader);
            parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
        }
        FileFormat::SdfV3000 => super::sdf_v3000::parse_sdf_v3000_file(path),
        FileFormat::Mol2 => super::mol2::parse_mol2_file(path),
        FileFormat::Xyz => super::xyz::parse_xyz_file(path),
    }
}

/// Parses multiple molecules from a file with automatic format detection.
///
/// This function reads the file, detects whether it is SDF V2000, V3000,
/// MOL2, or XYZ format, and uses the appropriate parser.
pub fn parse_auto_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    let content = std::fs::read_to_string(&path)?;
    match detect_format(&content) {
        FileFormat::SdfV2000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            let iter = SdfIterator::new(reader);
            iter.collect()
        }
        FileFormat::SdfV3000 => super::sdf_v3000::parse_sdf_v3000_file_multi(path),
        FileFormat::Mol2 => super::mol2::parse_mol2_file_multi(path),
        FileFormat::Xyz => super::xyz::parse_xyz_file_multi(path),
    }
}

/// A boxed iterator over molecules that works with any supported format.
pub type AutoIterator = Box<dyn Iterator<Item = Result<Molecule>>>;

/// Returns an iterator over molecules in a file with automatic format detection.
///
/// This function reads the file, detects whether it is SDF V2000, V3000,
/// MOL2, or XYZ format, and returns an appropriate iterator.
///
/// # Example
///
/// ```rust,ignore
/// use sdfrust::iter_auto_file;
///
/// // Iterate over any supported format
/// for result in iter_auto_file("large_database.mol2")? {
///     let mol = result?;
///     process(mol);
/// }
/// ```
pub fn iter_auto_file<P: AsRef<std::path::Path>>(path: P) -> Result<AutoIterator> {
    // Read enough to detect format (first few lines)
    let content = std::fs::read_to_string(&path)?;
    let format = detect_format(&content);

    match format {
        FileFormat::SdfV2000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            Ok(Box::new(SdfIterator::new(reader)))
        }
        FileFormat::SdfV3000 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            Ok(Box::new(super::sdf_v3000::SdfV3000Iterator::new(reader)))
        }
        FileFormat::Mol2 => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            Ok(Box::new(super::mol2::Mol2Iterator::new(reader)))
        }
        FileFormat::Xyz => {
            let cursor = std::io::Cursor::new(content);
            let reader = std::io::BufReader::new(cursor);
            Ok(Box::new(super::xyz::XyzIterator::new(reader)))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_MOL: &str = r#"methane
  test    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;

    #[test]
    fn test_parse_simple_molecule() {
        let mol = parse_sdf_string(SIMPLE_MOL).unwrap();

        assert_eq!(mol.name, "methane");
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 4);
        assert_eq!(mol.formula(), "CH4");

        // Check first atom (carbon)
        let carbon = &mol.atoms[0];
        assert_eq!(carbon.element, "C");
        assert_eq!(carbon.x, 0.0);
        assert_eq!(carbon.y, 0.0);
        assert_eq!(carbon.z, 0.0);

        // Check all bonds are single bonds
        for bond in &mol.bonds {
            assert_eq!(bond.order, BondOrder::Single);
        }
    }

    #[test]
    fn test_parse_with_properties() {
        let mol_with_props = r#"aspirin
  test    3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
M  END
> <MW>
180.16

> <SMILES>
CC(=O)OC1=CC=CC=C1C(=O)O

$$$$
"#;
        let mol = parse_sdf_string(mol_with_props).unwrap();

        assert_eq!(mol.name, "aspirin");
        assert_eq!(mol.get_property("MW"), Some("180.16"));
        assert_eq!(mol.get_property("SMILES"), Some("CC(=O)OC1=CC=CC=C1C(=O)O"));
    }

    #[test]
    fn test_multi_molecule_parsing() {
        let multi_mol = format!("{}{}", SIMPLE_MOL, SIMPLE_MOL);
        let mols = parse_sdf_string_multi(&multi_mol).unwrap();

        assert_eq!(mols.len(), 2);
        assert_eq!(mols[0].name, "methane");
        assert_eq!(mols[1].name, "methane");
    }
}
