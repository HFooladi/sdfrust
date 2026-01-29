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
        "ar" | "am" => BondOrder::Aromatic,
        "du" => BondOrder::Single, // Dummy bond
        "un" => BondOrder::Single, // Unknown, treat as single
        "nc" => BondOrder::Single, // Not connected, treat as single
        _ => BondOrder::Single,    // Default to single
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

        // Line 4: charge type (optional, skip)
        if !self.read_line()? {
            return Err(SdfError::MissingSection("MOLECULE charge type".to_string()));
        }

        // Skip remaining lines until next section
        // (status bits, comment, etc.)

        let mut atoms = Vec::with_capacity(num_atoms);
        let mut bonds = Vec::with_capacity(num_bonds);
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
                            atoms = self.parse_atom_section(num_atoms)?;
                        }
                        "BOND" => {
                            bonds = self.parse_bond_section(num_bonds, atoms.len())?;
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
    fn parse_atom_section(&mut self, expected_count: usize) -> Result<Vec<Atom>> {
        let mut atoms = Vec::with_capacity(expected_count);

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

            let atom = self.parse_atom_line(i)?;
            atoms.push(atom);
        }

        Ok(atoms)
    }

    /// Parses a single atom line.
    /// Format: atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    fn parse_atom_line(&self, index: usize) -> Result<Atom> {
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
        let element = atom_type.split('.').next().unwrap_or(atom_type).to_string();

        // Parse partial charge if present
        let partial_charge: f64 = if parts.len() > 8 {
            parts[8].parse().unwrap_or(0.0)
        } else {
            0.0
        };

        // Round partial charge to formal charge (rough approximation)
        let formal_charge = partial_charge.round() as i8;

        Ok(Atom {
            index,
            element,
            x,
            y,
            z,
            formal_charge,
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

    /// Parses the BOND section.
    fn parse_bond_section(
        &mut self,
        expected_count: usize,
        atom_count: usize,
    ) -> Result<Vec<Bond>> {
        let mut bonds = Vec::with_capacity(expected_count);

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

            let bond = self.parse_bond_line(atom_count)?;
            bonds.push(bond);
        }

        Ok(bonds)
    }

    /// Parses a single bond line.
    /// Format: bond_id origin_atom_id target_atom_id bond_type [status_bits]
    fn parse_bond_line(&self, atom_count: usize) -> Result<Bond> {
        let parts: Vec<&str> = self.current_line.split_whitespace().collect();

        if parts.len() < 4 {
            return Err(SdfError::Parse {
                line: self.line_number,
                message: format!("Bond line too short: {}", self.current_line),
            });
        }

        // parts[0] = bond_id (ignored)
        // parts[1] = origin_atom_id (1-based)
        // parts[2] = target_atom_id (1-based)
        // parts[3] = bond_type

        let atom1: usize = parts[1]
            .parse::<usize>()
            .map_err(|_| SdfError::Parse {
                line: self.line_number,
                message: "Invalid atom1 index".to_string(),
            })?
            .checked_sub(1)
            .ok_or(SdfError::InvalidAtomIndex {
                index: 0,
                atom_count,
            })?;

        let atom2: usize = parts[2]
            .parse::<usize>()
            .map_err(|_| SdfError::Parse {
                line: self.line_number,
                message: "Invalid atom2 index".to_string(),
            })?
            .checked_sub(1)
            .ok_or(SdfError::InvalidAtomIndex {
                index: 0,
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
                self.finished = true;
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
pub fn parse_mol2_file<P: AsRef<std::path::Path>>(path: P) -> Result<Molecule> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let mut parser = Mol2Parser::new(reader);

    parser.parse_molecule()?.ok_or(SdfError::EmptyFile)
}

/// Parses all molecules from a MOL2 file.
pub fn parse_mol2_file_multi<P: AsRef<std::path::Path>>(path: P) -> Result<Vec<Molecule>> {
    let file = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(file);
    let iter = Mol2Iterator::new(reader);

    iter.collect()
}

/// Returns an iterator over molecules in a MOL2 file.
pub fn iter_mol2_file<P: AsRef<std::path::Path>>(
    path: P,
) -> Result<Mol2Iterator<std::io::BufReader<std::fs::File>>> {
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
