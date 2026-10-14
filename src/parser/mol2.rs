//! TRIPOS MOL2 format parser.
//!
//! MOL2 is a common format for storing molecular structures with atom types
//! and partial charges. It uses a record-based format with sections marked
//! by `@<TRIPOS>SECTION_NAME`.
//!
//! ## Supported Sections
//! - `@<TRIPOS>MOLECULE` - Molecule name and counts
//! - `@<TRIPOS>ATOM` - Atom coordinates, types, and charges
//! - `@<TRIPOS>BOND` - Bond connectivity and types

use std::collections::HashMap;
use std::io::BufRead;

use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::error::{Result, SdfError};
use crate::molecule::Molecule;

/// MOL2 bond type to BondOrder mapping
fn mol2_bond_type_to_order(bond_type: &str) -> BondOrder {
    match bond_type {
        "1" | "1n" => BondOrder::Single,
        "2" => BondOrder::Double,
        "3" => BondOrder::Triple,
        "ar" => BondOrder::Aromatic,
        "am" => BondOrder::Single, // Amide: a single bond with partial double character
        "du" => BondOrder::Single, // Dummy bond
        "un" => BondOrder::Single, // Unknown, treat as single
        "nc" => BondOrder::Single, // Not connected, treat as single
        _ => BondOrder::Single,    // Default to single
    }
}

/// Extracts and normalizes the element from a SYBYL atom type
/// (e.g. "C.ar" → "C", "CL" → "Cl"). Pseudo-atoms (Du, LP, Any, ...) are kept as-is.
fn element_from_atom_type(atom_type: &str) -> String {
    let symbol = atom_type.split('.').next().unwrap_or(atom_type);
    match symbol {
        "Du" | "LP" | "Lp" | "Any" | "Hal" | "Het" | "Hev" => symbol.to_string(),
        _ => {
            let mut chars = symbol.chars();
            match chars.next() {
                Some(first) => {
                    first.to_ascii_uppercase().to_string() + &chars.as_str().to_ascii_lowercase()
                }
                None => String::new(),
            }
        }
    }
}

/// TRIPOS MOL2 format parser.
pub struct Mol2Parser<R> {
    reader: R,
    line_number: usize,
    current_line: String,
    peeked: bool,
}

impl<R: BufRead> Mol2Parser<R> {
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

    /// Peeks at the current line without consuming it.
    fn peek_line(&mut self) -> Result<Option<&str>> {
        if !self.peeked {
            if !self.read_line()? {
                return Ok(None);
            }
            self.peeked = true;
        }
        Ok(Some(&self.current_line))
    }

    /// Skips lines until a section header is found.
    fn skip_to_section(&mut self) -> Result<Option<String>> {
        loop {
            if !self.read_line()? {
                return Ok(None);
            }
            if self.current_line.starts_with("@<TRIPOS>") {
                let section = self.current_line[9..].trim().to_uppercase();
                return Ok(Some(section));
            }
        }
    }

    /// Parses a single molecule from the input.
    /// Returns None if end of file is reached.
    pub fn parse_molecule(&mut self) -> Result<Option<Molecule>> {
        // Find @<TRIPOS>MOLECULE section
        loop {
            match self.skip_to_section()? {
                Some(section) if section == "MOLECULE" => break,
                Some(_) => continue,
                None => return Ok(None),
            }
        }

        // Parse MOLECULE section
        // Line 1: molecule name
        if !self.read_line()? {
            return Err(SdfError::MissingSection("MOLECULE name".to_string()));
        }
        let name = self.current_line.trim().to_string();

        // Line 2: num_atoms num_bonds [num_subst num_feat num_sets]
        if !self.read_line()? {
            return Err(SdfError::MissingSection("MOLECULE counts".to_string()));
        }
        let counts: Vec<&str> = self.current_line.split_whitespace().collect();
        if counts.len() < 2 {
            return Err(SdfError::InvalidCountsLine(self.current_line.clone()));
        }

        let num_atoms: usize = counts[0].parse().map_err(|_| {
            SdfError::InvalidCountsLine(format!("Invalid atom count: {}", counts[0]))
        })?;
        let num_bonds: usize = counts[1].parse().map_err(|_| {
            SdfError::InvalidCountsLine(format!("Invalid bond count: {}", counts[1]))
        })?;

        // Line 3: molecule type (optional, skip)
        if !self.read_line()? {
            return Err(SdfError::MissingSection("MOLECULE type".to_string()));
        }

        // Line 4: charge type
        if !self.read_line()? {
            return Err(SdfError::MissingSection("MOLECULE charge type".to_string()));
        }
        let has_charges = !self.current_line.trim().eq_ignore_ascii_case("NO_CHARGES");

        // Skip remaining lines until next section
        // (status bits, comment, etc.)

        let mut atoms = Vec::new();
        let mut atom_ids: HashMap<String, usize> = HashMap::new();
        let mut bonds = Vec::new();
        let properties = HashMap::new();

        // Parse sections until we hit another MOLECULE or EOF
        loop {
            match self.peek_line()? {
                None => break,
                Some(line) if line.starts_with("@<TRIPOS>MOLECULE") => break,
                Some(line) if line.starts_with("@<TRIPOS>") => {
                    let section = line[9..].trim().to_uppercase();
                    self.read_line()?; // consume the section header

                    match section.as_str() {
                        "ATOM" => {
                            (atoms, atom_ids) = self.parse_atom_section(num_atoms, has_charges)?;
                        }
                        "BOND" => {
                            bonds = self.parse_bond_section(num_bonds, &atom_ids)?;
                        }
                        "SUBSTRUCTURE" | "COMMENT" | "HEADTAIL" | "SET" | "DICT" | "EXTENSION" => {
                            // Skip these sections
                            self.skip_section()?;
                        }
                        _ => {
                            // Unknown section, skip
                            self.skip_section()?;
                        }
                    }
                }
                _ => {
                    self.read_line()?; // consume non-section line
                }
            }
        }

        // Validate counts
        if atoms.len() != num_atoms {
            return Err(SdfError::AtomCountMismatch {
                expected: num_atoms,
                found: atoms.len(),
            });
        }

        Ok(Some(Molecule {
            name,
            program_line: None,
            comment: None,
            atoms,
            bonds,
            properties,
            format_version: crate::molecule::SdfFormat::V2000,
            stereogroups: Vec::new(),
            sgroups: Vec::new(),
            collections: Vec::new(),
        }))
    }

    /// Parses the ATOM section.
    fn parse_atom_section(
        &mut self,
        expected_count: usize,
        has_charges: bool,
    ) -> Result<(Vec<Atom>, HashMap<String, usize>)> {
        let prealloc = expected_count.min(super::sdf::MAX_PREALLOC);
        let mut atoms = Vec::with_capacity(prealloc);
        let mut atom_ids = HashMap::with_capacity(prealloc);

        for i in 0..expected_count {
            if !self.read_line()? {
                return Err(SdfError::AtomCountMismatch {
                    expected: expected_count,
                    found: i,
                });
            }

            // Check if we hit another section
            if self.current_line.starts_with("@<TRIPOS>") {
                self.peeked = true;
                return Err(SdfError::AtomCountMismatch {
                    expected: expected_count,
                    found: i,
                });
            }

            let atom = self.parse_atom_line(i, has_charges)?;
            let id = self
                .current_line
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_string();
            if atom_ids.insert(id.clone(), i).is_some() {
                return Err(SdfError::Parse {
                    line: self.line_number,
                    message: format!("Duplicate MOL2 atom id: {}", id),
                });
            }
            atoms.push(atom);
        }

        Ok((atoms, atom_ids))
    }

    /// Parses a single atom line.
    /// Format: atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    fn parse_atom_line(&self, index: usize, has_charges: bool) -> Result<Atom> {
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();

        if parts.len() < 6 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Atom line too short: {}", self.current_line),
            });
        }

        // atom_id is 1-based, we use 0-based index
        // parts[0] = atom_id (ignored, we use our own index)
        // parts[1] = atom_name
        // parts[2] = x
        // parts[3] = y
        // parts[4] = z
        // parts[5] = atom_type (e.g., "C.3", "N.ar", "O.2")
        // parts[6] = subst_id (optional)
        // parts[7] = subst_name (optional)
        // parts[8] = charge (optional)

        let x: f64 = parts[2]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[2].to_string()))?;
        let y: f64 = parts[3]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[3].to_string()))?;
        let z: f64 = parts[4]
            .parse()
            .map_err(|_| SdfError::InvalidCoordinate(parts[4].to_string()))?;

        // Extract element from atom_type (e.g., "C.3" -> "C", "N.ar" -> "N")
        let atom_type = parts[5];
        let element = element_from_atom_type(atom_type);

        // Partial charge (kept separate: it is not a formal charge)
        let partial_charge: Option<f64> = if has_charges && parts.len() > 8 {
            Some(parts[8].parse().map_err(|_| SdfError::Parse {
                line: self.line_number,
                message: format!("Invalid partial charge: {}", parts[8]),
            })?)
        } else {
            None
        };

        // MOL2 has no formal-charge column; the only unambiguous case is the
        // positively charged sp3 nitrogen type.
        let formal_charge = if atom_type == "N.4" { 1 } else { 0 };

        Ok(Atom {
            index,
            element,
            x,
            y,
            z,
            formal_charge,
            mass_difference: 0,
            isotope: None,
            partial_charge,
            no_implicit: true,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
            v3000_id: None,
            atom_atom_mapping: None,
            rgroup_label: None,
            radical: None,
        })
    }

    /// Parses the BOND section.
    fn parse_bond_section(
        &mut self,
        expected_count: usize,
        atom_ids: &HashMap<String, usize>,
    ) -> Result<Vec<Bond>> {
        let mut bonds = Vec::with_capacity(expected_count.min(super::sdf::MAX_PREALLOC));

        for i in 0..expected_count {
            if !self.read_line()? {
                return Err(SdfError::BondCountMismatch {
                    expected: expected_count,
                    found: i,
                });
            }

            // Check if we hit another section
            if self.current_line.starts_with("@<TRIPOS>") {
                self.peeked = true;
                return Err(SdfError::BondCountMismatch {
                    expected: expected_count,
                    found: i,
                });
            }

            let bond = self.parse_bond_line(atom_ids)?;
            bonds.push(bond);
        }

        Ok(bonds)
    }

    /// Parses a single bond line.
    /// Format: bond_id origin_atom_id target_atom_id bond_type [status_bits]
    fn parse_bond_line(&self, atom_ids: &HashMap<String, usize>) -> Result<Bond> {
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();

        if parts.len() < 4 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Bond line too short: {}", self.current_line),
            });
        }

        // parts[0] = bond_id (ignored)
        // parts[1] = origin_atom_id, parts[2] = target_atom_id (atom ids, which
        // need not be contiguous), parts[3] = bond_type
        let resolve = |id: &str| {
            atom_ids
                .get(id)
                .copied()
                .ok_or_else(|| SdfError::InvalidAtomIndex {
                    index: id.parse().unwrap_or(0),
                    atom_count: atom_ids.len(),
                })
        };
        let atom1 = resolve(parts[1])?;
        let atom2 = resolve(parts[2])?;

        let bond_type = parts[3].to_lowercase();
        let order = mol2_bond_type_to_order(&bond_type);

        Ok(Bond {
            atom1,
            atom2,
            order,
            stereo: BondStereo::None,
            topology: None,
            v3000_id: None,
            reacting_center: None,
        })
    }

    /// Skips lines until the next section header.
    fn skip_section(&mut self) -> Result<()> {
        loop {
            match self.peek_line()? {
                None => return Ok(()),
                Some(line) if line.starts_with("@<TRIPOS>") => return Ok(()),
                _ => {
                    self.read_line()?;
                }
            }
        }
    }
}

/// Iterator over molecules in a MOL2 file.
pub struct Mol2Iterator<R> {
    parser: Mol2Parser<R>,
    finished: bool,
}

impl<R: BufRead> Mol2Iterator<R> {
    /// Creates a new iterator from a buffered reader.
    pub fn new(reader: R) -> Self {
        Self {
            parser: Mol2Parser::new(reader),
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for Mol2Iterator<R> {
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
                // The next call skips ahead to the next @<TRIPOS>MOLECULE, so only
                // I/O errors end the iteration.
                if matches!(e, SdfError::Io(_)) {
                    self.finished = true;
                }
                Some(Err(e))
            }
        }
    }
}

/// Parses a single molecule from a MOL2 string.
pub fn parse_mol2_string(content: &str) -> Result<Molecule> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let mut parser = Mol2Parser::new(reader);

    parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
}

/// Parses all molecules from a MOL2 string.
pub fn parse_mol2_string_multi(content: &str) -> Result<Vec<Molecule>> {
    let cursor = std::io::Cursor::new(content);
    let reader = std::io::BufReader::new(cursor);
    let iter = Mol2Iterator::new(reader);

    iter.collect()
}

/// Parses a single molecule from a MOL2 file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed.
pub fn parse_mol2_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    #[cfg(feature = "gzip")]
    {
        let reader = super::compression::open_maybe_gz(&path)?;
        let mut parser = Mol2Parser::new(reader);
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
        let mut parser = Mol2Parser::new(reader);
        parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
    }
}

/// Parses all molecules from a MOL2 file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed.
pub fn parse_mol2_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    #[cfg(feature = "gzip")]
    {
        let reader = super::compression::open_maybe_gz(&path)?;
        let iter = Mol2Iterator::new(reader);
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
        let iter = Mol2Iterator::new(reader);
        iter.collect()
    }
}

/// Returns an iterator over molecules in a MOL2 file.
///
/// When the `gzip` feature is enabled, files ending in `.gz` are automatically
/// decompressed. Note that the return type differs based on the feature flag.
#[cfg(feature = "gzip")]
pub fn iter_mol2_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<Mol2Iterator<super::compression::MaybeGzReader>> {
    let reader = super::compression::open_maybe_gz(&path)?;
    Ok(Mol2Iterator::new(reader))
}

/// Returns an iterator over molecules in a MOL2 file.
#[cfg(not(feature = "gzip"))]
pub fn iter_mol2_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<Mol2Iterator<std::io::BufReader<std::fs::File>>> {
    if path
        .as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
    {
        return Err(SdfError::GzipNotEnabled);
    }
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    Ok(Mol2Iterator::new(reader))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_MOL2: &str = r#"@<TRIPOS>MOLECULE
methane
 5 4 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 H1          0.6289    0.6289    0.6289 H         1 MOL       0.0000
      3 H2         -0.6289   -0.6289    0.6289 H         1 MOL       0.0000
      4 H3         -0.6289    0.6289   -0.6289 H         1 MOL       0.0000
      5 H4          0.6289   -0.6289   -0.6289 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
"#;

    #[test]
    fn test_parse_simple_mol2() {
        let mol = parse_mol2_string(SIMPLE_MOL2).unwrap();

        assert_eq!(mol.name, "methane");
        assert_eq!(mol.atom_count(), 5);
        assert_eq!(mol.bond_count(), 4);
        assert_eq!(mol.formula(), "CH4");
    }

    #[test]
    fn test_parse_mol2_atoms() {
        let mol = parse_mol2_string(SIMPLE_MOL2).unwrap();

        // Check carbon
        assert_eq!(mol.atoms[0].element, "C");
        assert_eq!(mol.atoms[0].x, 0.0);

        // Check hydrogen
        assert_eq!(mol.atoms[1].element, "H");
    }

    #[test]
    fn test_parse_mol2_bonds() {
        let mol = parse_mol2_string(SIMPLE_MOL2).unwrap();

        // All bonds should be single
        for bond in &mol.bonds {
            assert_eq!(bond.order, BondOrder::Single);
        }

        // First bond connects C to H
        assert_eq!(mol.bonds[0].atom1, 0);
        assert_eq!(mol.bonds[0].atom2, 1);
    }

    #[test]
    fn test_parse_mol2_aromatic() {
        let benzene_mol2 = r#"@<TRIPOS>MOLECULE
benzene
 6 6 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          1.2124    0.7000    0.0000 C.ar      1 MOL       0.0000
      2 C2          1.2124   -0.7000    0.0000 C.ar      1 MOL       0.0000
      3 C3          0.0000   -1.4000    0.0000 C.ar      1 MOL       0.0000
      4 C4         -1.2124   -0.7000    0.0000 C.ar      1 MOL       0.0000
      5 C5         -1.2124    0.7000    0.0000 C.ar      1 MOL       0.0000
      6 C6          0.0000    1.4000    0.0000 C.ar      1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 ar
     2     2     3 ar
     3     3     4 ar
     4     4     5 ar
     5     5     6 ar
     6     6     1 ar
"#;

        let mol = parse_mol2_string(benzene_mol2).unwrap();

        assert_eq!(mol.name, "benzene");
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);

        // All bonds should be aromatic
        for bond in &mol.bonds {
            assert_eq!(bond.order, BondOrder::Aromatic);
        }
    }

    #[test]
    fn test_parse_mol2_with_charges() {
        let charged_mol2 = r#"@<TRIPOS>MOLECULE
ammonium
 5 4 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 N1          0.0000    0.0000    0.0000 N.4       1 MOL       1.0000
      2 H1          0.6289    0.6289    0.6289 H         1 MOL       0.0000
      3 H2         -0.6289   -0.6289    0.6289 H         1 MOL       0.0000
      4 H3         -0.6289    0.6289   -0.6289 H         1 MOL       0.0000
      5 H4          0.6289   -0.6289   -0.6289 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
"#;

        let mol = parse_mol2_string(charged_mol2).unwrap();

        assert_eq!(mol.name, "ammonium");
        assert_eq!(mol.atoms[0].formal_charge, 1);
        assert_eq!(mol.total_charge(), 1);
    }

    #[test]
    fn test_parse_multi_mol2() {
        let multi = format!("{}{}", SIMPLE_MOL2, SIMPLE_MOL2);
        let mols = parse_mol2_string_multi(&multi).unwrap();

        assert_eq!(mols.len(), 2);
        assert_eq!(mols[0].name, "methane");
        assert_eq!(mols[1].name, "methane");
    }
}
