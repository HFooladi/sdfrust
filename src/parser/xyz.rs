//! XYZ format parser.
//!
//! XYZ is a simple format for molecular coordinates containing only atoms
//! (no bond information). It is commonly used in computational chemistry.
//!
//! ## Format
//!
//! ```text
//! 3                        <- atom count (integer)
//! water molecule           <- comment/title (molecule name)
//! O  0.000000  0.000000  0.117300    <- element x y z
//! H  0.756950  0.000000 -0.469200
//! H -0.756950  0.000000 -0.469200
//! ```
//!
//! Multi-molecule files concatenate blocks (no delimiter).

use std::collections::HashMap;
use std::io::BufRead;

use crate::atom::Atom;
use crate::error::{Result, SdfError};
use crate::molecule::Molecule;

/// Maps atomic numbers to element symbols.
fn atomic_number_to_symbol(num: u8) -> Option<&'static str> {
    match num {
        1 => Some("H"),
        2 => Some("He"),
        3 => Some("Li"),
        4 => Some("Be"),
        5 => Some("B"),
        6 => Some("C"),
        7 => Some("N"),
        8 => Some("O"),
        9 => Some("F"),
        10 => Some("Ne"),
        11 => Some("Na"),
        12 => Some("Mg"),
        13 => Some("Al"),
        14 => Some("Si"),
        15 => Some("P"),
        16 => Some("S"),
        17 => Some("Cl"),
        18 => Some("Ar"),
        19 => Some("K"),
        20 => Some("Ca"),
        26 => Some("Fe"),
        29 => Some("Cu"),
        30 => Some("Zn"),
        35 => Some("Br"),
        53 => Some("I"),
        _ => None,
    }
}

/// Normalizes an element symbol to proper case (first letter uppercase, rest lowercase).
fn normalize_element(element: &str) -> String {
    let element = element.trim();
    if element.is_empty() {
        return String::new();
    }

    let mut chars = element.chars();
    let first = chars.next().unwrap().to_uppercase().to_string();
    let rest: String = chars.collect::<String>().to_lowercase();
    first + &rest
}

/// XYZ format parser.
pub struct XyzParser<R> {
    reader: R,
    line_number: usize,
    current_line: String,
    peeked: bool,
}

impl<R: BufRead> XyzParser<R> {
    /// Creates a new parser from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
            current_line: String::new(),
            peeked: false,
        }
    }

    /// Reads the next line from the input.
    fn read_line(&mut self) -> Result<bool> {
        if self.peeked {
            self.peeked = false;
            return Ok(!self.current_line.is_empty());
        }

        self.current_line.clear();
        let bytes_read = self.reader.read_line(&mut self.current_line)?;
        if bytes_read > 0 {
            self.line_number += 1;
            // Remove trailing newline
            if self.current_line.ends_with('\n') {
                self.current_line.pop();
                if self.current_line.ends_with('\r') {
                    self.current_line.pop();
                }
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Skips blank lines and returns true if a non-blank line was found.
    fn skip_blank_lines(&mut self) -> Result<bool> {
        loop {
            if !self.read_line()? {
                return Ok(false);
            }
            if !self.current_line.trim().is_empty() {
                self.peeked = true;
                return Ok(true);
            }
        }
    }

    /// Parses the atom count line (first line of XYZ block).
    fn parse_atom_count_line(&self) -> Result<usize> {
        let trimmed = self.current_line.trim();
        trimmed
            .parse::<usize>()
            .map_err(|_| SdfError::InvalidCountsLine(format!("Invalid atom count: {}", trimmed)))
    }

    /// Parses a single atom line.
    /// Format: element x y z [additional columns ignored]
    fn parse_atom_line(&self, index: usize) -> Result<Atom> {
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();

        if parts.len() < 4 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Atom line too short: {}", self.current_line),
            });
        }

        // First field is element (could be symbol like "C" or atomic number like "6")
        let element_str = parts[0];
        let element = if let Ok(atomic_num) = element_str.parse::<u8>() {
            // It's an atomic number
            atomic_number_to_symbol(atomic_num)
                .map(|s| s.to_string())
                .ok_or_else(|| SdfError::Parse {
                    line: self.line_number,
                    message: format!("Unknown atomic number: {}", atomic_num),
                })?
        } else {
            // It's an element symbol - normalize case
            normalize_element(element_str)
        };

        let x: f64 = parts[1]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[1].to_string()))?;
        let y: f64 = parts[2]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[2].to_string()))?;
        let z: f64 = parts[3]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[3].to_string()))?;

        Ok(Atom {
            index,
            element,
            x,
            y,
            z,
            formal_charge: 0,
            mass_difference: 0,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
            v3000_id: None,
            atom_atom_mapping: None,
            rgroup_label: None,
            radical: None,
        })
    }

    /// Parses a single molecule from the input.
    /// Returns None if end of file is reached.
    pub fn parse_molecule(&mut self) -> Result<Option<Molecule>> {
        // Skip any leading blank lines
        if !self.skip_blank_lines()? {
            return Ok(None);
        }

        // Line 1: Atom count
        if !self.read_line()? {
            return Ok(None);
        }
        let num_atoms = self.parse_atom_count_line()?;

        // Line 2: Comment/title (molecule name)
        if !self.read_line()? {
            return Err(SdfError::MissingSection("XYZ comment line".to_string()));
        }
        let name = self.current_line.trim().to_string();

        // Parse atoms
        let mut atoms = Vec::with_capacity(num_atoms);
        for i in 0..num_atoms {
            if !self.read_line()? {
                return Err(SdfError::AtomCountMismatch {
                    expected: num_atoms,
                    found: i,
                });
            }
            let atom = self.parse_atom_line(i)?;
            atoms.push(atom);
        }

        Ok(Some(Molecule {
            name,
            program_line: None,
            comment: None,
            atoms,
            bonds: Vec::new(), // XYZ has no bond information
            properties: HashMap::new(),
            format_version: crate::molecule::SdfFormat::V2000,
            stereogroups: Vec::new(),
            sgroups: Vec::new(),
            collections: Vec::new(),
        }))
    }
}

/// Iterator over molecules in an XYZ file.
pub struct XyzIterator<R> {
    parser: XyzParser<R>,
    finished: bool,
}

impl<R: BufRead> XyzIterator<R> {
    /// Creates a new iterator from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            parser: XyzParser::new(reader),
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for XyzIterator<R> {
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

/// Parses a single molecule from an XYZ string.
pub fn parse_xyz_string(content: &str) -> Result<Molecule> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let mut parser = XyzParser::new(reader);

    parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
}

/// Parses all molecules from an XYZ string.
pub fn parse_xyz_string_multi(content: &str) -> Result<Vec<Molecule>> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let iter = XyzIterator::new(reader);

    iter.collect()
}

/// Parses a single molecule from an XYZ file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed.
pub fn parse_xyz_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    #[cfg(feature = "gzip")]
    {
        let reader = super::compression::open_maybe_gz(&path)?;
        let mut parser = XyzParser::new(reader);
        parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
    }

    #[cfg(not(feature = "gzip"))]
    {
        if path
            .as_ref()
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            return Err(SdfError::GzipNotEnabled);
        }
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        let mut parser = XyzParser::new(reader);
        parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
    }
}

/// Parses all molecules from an XYZ file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed.
pub fn parse_xyz_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    #[cfg(feature = "gzip")]
    {
        let reader = super::compression::open_maybe_gz(&path)?;
        let iter = XyzIterator::new(reader);
        iter.collect()
    }

    #[cfg(not(feature = "gzip"))]
    {
        if path
            .as_ref()
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            return Err(SdfError::GzipNotEnabled);
        }
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        let iter = XyzIterator::new(reader);
        iter.collect()
    }
}

/// Returns an iterator over molecules in an XYZ file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed. Note that the return type differs based on the feature flag.
#[cfg(feature = "gzip")]
pub fn iter_xyz_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<XyzIterator<super::compression::MaybeGzReader>> {
    let reader = super::compression::open_maybe_gz(&path)?;
    Ok(XyzIterator::new(reader))
}

/// Returns an iterator over molecules in an XYZ file.
#[cfg(not(feature = "gzip"))]
pub fn iter_xyz_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<XyzIterator<std::io::BufReader<std::fs::File>>> {
    if path
        .as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
    {
        return Err(SdfError::GzipNotEnabled);
    }
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    Ok(XyzIterator::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_XYZ: &str = r#"3
water molecule
O  0.000000  0.000000  0.117300
H  0.756950  0.000000 -0.469200
H -0.756950  0.000000 -0.469200
"#;

    const METHANE_XYZ: &str = r#"5
methane
C   0.000000   0.000000   0.000000
H   0.628900   0.628900   0.628900
H  -0.628900  -0.628900   0.628900
H  -0.628900   0.628900  -0.628900
H   0.628900  -0.628900  -0.628900
"#;

    #[test]
    fn test_parse_water() {
        let mol = parse_xyz_string(WATER_XYZ).unwrap();

        assert_eq!(mol.name, "water molecule");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 0); // XYZ has no bonds
        assert_eq!(mol.formula(), "H2O");

        // Check oxygen
        assert_eq!(mol.atoms[0].element, "O");
        assert!((mol.atoms[0].x - 0.0).abs() < 1e-6);
        assert!((mol.atoms[0].y - 0.0).abs() < 1e-6);
        assert!((mol.atoms[0].z - 0.1173).abs() < 1e-6);

        // Check first hydrogen
        assert_eq!(mol.atoms[1].element, "H");
        assert!((mol.atoms[1].x - 0.75695).abs() < 1e-6);
    }

    #[test]
    fn test_parse_methane() {
        let mol = parse_xyz_string(METHANE_XYZ).unwrap();

        assert_eq!(mol.name, "methane");
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.formula(), "CH4");
    }

    #[test]
    fn test_parse_multi() {
        let multi = format!("{}{}", WATER_XYZ, METHANE_XYZ);
        let mols = parse_xyz_string_multi(&multi).unwrap();

        assert_eq!(mols.len(), 2);
        assert_eq!(mols[0].name, "water molecule");
        assert_eq!(mols[1].name, "methane");
    }

    #[test]
    fn test_parse_atomic_numbers() {
        let xyz = r#"3
atomic number test
8  0.0  0.0  0.0
1  1.0  0.0  0.0
1 -1.0  0.0  0.0
"#;
        let mol = parse_xyz_string(xyz).unwrap();

        assert_eq!(mol.atoms[0].element, "O");
        assert_eq!(mol.atoms[1].element, "H");
        assert_eq!(mol.atoms[2].element, "H");
    }

    #[test]
    fn test_normalize_element_case() {
        let xyz = r#"3
case test
o  0.0  0.0  0.0
h  1.0  0.0  0.0
CA -1.0  0.0  0.0
"#;
        let mol = parse_xyz_string(xyz).unwrap();

        assert_eq!(mol.atoms[0].element, "O");
        assert_eq!(mol.atoms[1].element, "H");
        assert_eq!(mol.atoms[2].element, "Ca");
    }

    #[test]
    fn test_extra_columns_ignored() {
        let xyz = r#"2
extra columns test
C  0.0  0.0  0.0  0.5  extra_data
H  1.0  0.0  0.0  -0.2
"#;
        let mol = parse_xyz_string(xyz).unwrap();

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.atoms[0].element, "C");
        assert_eq!(mol.atoms[1].element, "H");
    }

    #[test]
    fn test_blank_lines_between_molecules() {
        let xyz = format!("{}\n\n{}", WATER_XYZ, METHANE_XYZ);
        let mols = parse_xyz_string_multi(&xyz).unwrap();

        assert_eq!(mols.len(), 2);
    }

    #[test]
    fn test_empty_file() {
        let result = parse_xyz_string("");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), SdfError::EmptyFile));
    }

    #[test]
    fn test_invalid_atom_count() {
        let xyz = r#"abc
test
C  0.0  0.0  0.0
"#;
        let result = parse_xyz_string(xyz);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SdfError::InvalidCountsLine(_)
        ));
    }

    #[test]
    fn test_fewer_atoms_than_declared() {
        let xyz = r#"5
missing atoms
C  0.0  0.0  0.0
H  1.0  0.0  0.0
"#;
        let result = parse_xyz_string(xyz);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SdfError::AtomCountMismatch {
                expected: 5,
                found: 2
            }
        ));
    }

    #[test]
    fn test_invalid_coordinate() {
        let xyz = r#"1
bad coords
C  abc  0.0  0.0
"#;
        let result = parse_xyz_string(xyz);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SdfError::InvalidCoordinate(_)
        ));
    }

    #[test]
    fn test_iterator() {
        let multi = format!("{}{}{}", WATER_XYZ, METHANE_XYZ, WATER_XYZ);
        let cursor = std::io::Cursor::new(multi);
        let reader = std::io::BufReader::new(cursor);
        let iter = XyzIterator::new(reader);

        let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
        assert_eq!(mols.len(), 3);
        assert_eq!(mols[0].name, "water molecule");
        assert_eq!(mols[1].name, "methane");
        assert_eq!(mols[2].name, "water molecule");
    }
}
