//! TRIPOS MOL2 format writer.
//!
//! This module provides functions to write molecules to the MOL2 format.
//! The MOL2 format uses sections marked by `@<TRIPOS>SECTION_NAME`.

use std::io::Write;
use std::path::Path;

use crate::bond::BondOrder;
use crate::error::Result;
use crate::molecule::Molecule;

/// Converts a BondOrder to MOL2 bond type string.
fn order_to_mol2_bond_type(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "1",
        BondOrder::Double => "2",
        BondOrder::Triple => "3",
        BondOrder::Aromatic => "ar",
        BondOrder::SingleOrDouble => "1",
        BondOrder::SingleOrAromatic => "ar",
        BondOrder::DoubleOrAromatic => "ar",
        BondOrder::Any => "1",
        BondOrder::Coordination => "1",
        BondOrder::Hydrogen => "1",
    }
}

/// Infers a simple SYBYL atom type from element and bond information.
///
/// Since the parser discards SYBYL subtypes, we infer from:
/// - Element symbol
/// - Presence of aromatic bonds
/// - Bond count for hybridization hints
fn infer_sybyl_type(mol: &Molecule, atom_index: usize) -> String {
    let atom = &mol.atoms[atom_index];
    let element = &atom.element;

    // Check if atom has aromatic bonds
    let has_aromatic = mol
        .bonds
        .iter()
        .any(|b| (b.atom1 == atom_index || b.atom2 == atom_index) && b.is_aromatic());

    // Count bonds for hybridization hints
    let bond_count = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_index || b.atom2 == atom_index)
        .count();

    // Calculate sum of bond orders for this atom
    let bond_order_sum: f64 = mol
        .bonds
        .iter()
        .filter(|b| b.atom1 == atom_index || b.atom2 == atom_index)
        .map(|b| b.order.order())
        .sum();

    match element.as_str() {
        "C" => {
            if has_aromatic {
                "C.ar".to_string()
            } else if bond_order_sum > 3.5 {
                "C.1".to_string() // sp hybridization (triple bond)
            } else if bond_order_sum > 2.5 || bond_count == 3 {
                "C.2".to_string() // sp2 hybridization
            } else {
                "C.3".to_string() // sp3 hybridization
            }
        }
        "N" => {
            if has_aromatic {
                "N.ar".to_string()
            } else if bond_count == 4 || atom.formal_charge > 0 {
                "N.4".to_string() // quaternary nitrogen
            } else if bond_order_sum > 2.5 {
                "N.1".to_string() // sp
            } else if bond_order_sum > 1.5 || bond_count == 3 {
                "N.pl3".to_string() // planar sp2
            } else {
                "N.3".to_string() // sp3
            }
        }
        "O" => {
            if has_aromatic {
                "O.ar".to_string()
            } else if atom.formal_charge < 0 {
                "O.co2".to_string() // carboxylate oxygen
            } else if bond_order_sum > 1.5 {
                "O.2".to_string() // sp2 (carbonyl)
            } else {
                "O.3".to_string() // sp3 (hydroxyl, ether)
            }
        }
        "S" => {
            if has_aromatic {
                "S.ar".to_string()
            } else if bond_order_sum > 1.5 {
                "S.2".to_string()
            } else {
                "S.3".to_string()
            }
        }
        "P" => "P.3".to_string(),
        "H" => "H".to_string(),
        "F" => "F".to_string(),
        "Cl" => "Cl".to_string(),
        "Br" => "Br".to_string(),
        "I" => "I".to_string(),
        // Default: just use the element symbol
        _ => element.clone(),
    }
}

/// Writes a molecule to MOL2 format.
pub fn write_mol2<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
    // @<TRIPOS>MOLECULE section
    writeln!(writer, "@<TRIPOS>MOLECULE")?;
    writeln!(
        writer,
        "{}",
        if molecule.name.is_empty() {
            "unnamed"
        } else {
            &molecule.name
        }
    )?;
    writeln!(
        writer,
        " {} {} 0 0 0",
        molecule.atom_count(),
        molecule.bond_count()
    )?;
    writeln!(writer, "SMALL")?;
    writeln!(writer, "NO_CHARGES")?;
    writeln!(writer)?;

    // @<TRIPOS>ATOM section
    writeln!(writer, "@<TRIPOS>ATOM")?;
    for (i, atom) in molecule.atoms.iter().enumerate() {
        let atom_name = format!("{}{}", atom.element, i + 1);
        let atom_type = infer_sybyl_type(molecule, i);
        writeln!(
            writer,
            "{:>7} {:<4} {:>10.4} {:>10.4} {:>10.4} {:<6} {:>3} {:<4} {:>8.4}",
            i + 1,                     // atom_id (1-based)
            atom_name,                 // atom_name
            atom.x,                    // x coordinate
            atom.y,                    // y coordinate
            atom.z,                    // z coordinate
            atom_type,                 // SYBYL atom type
            1,                         // subst_id
            "MOL",                     // subst_name
            atom.formal_charge as f64  // charge (as float)
        )?;
    }

    // @<TRIPOS>BOND section
    writeln!(writer, "@<TRIPOS>BOND")?;
    for (i, bond) in molecule.bonds.iter().enumerate() {
        let bond_type = order_to_mol2_bond_type(bond.order);
        writeln!(
            writer,
            "{:>6} {:>5} {:>5} {}",
            i + 1,          // bond_id (1-based)
            bond.atom1 + 1, // origin_atom_id (1-based)
            bond.atom2 + 1, // target_atom_id (1-based)
            bond_type       // bond_type
        )?;
    }

    Ok(())
}

/// Writes a molecule to a MOL2 string.
pub fn write_mol2_string(molecule: &Molecule) -> Result<String> {
    let mut buffer = Vec::new();
    write_mol2(&mut buffer, molecule)?;
    Ok(String::from_utf8_lossy(&buffer).to_string())
}

/// Writes multiple molecules to MOL2 format.
pub fn write_mol2_multi<W: Write>(writer: &mut W, molecules: &[Molecule]) -> Result<()> {
    for mol in molecules {
        write_mol2(writer, mol)?;
    }
    Ok(())
}

/// Writes a molecule to a MOL2 file.
pub fn write_mol2_file<P: AsRef<Path>>(path: P, molecule: &Molecule) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_mol2(&mut writer, molecule)
}

/// Writes multiple molecules to a MOL2 file.
pub fn write_mol2_file_multi<P: AsRef<Path>>(path: P, molecules: &[Molecule]) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_mol2_multi(&mut writer, molecules)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;
    use crate::parser::parse_mol2_string;

    #[test]
    fn test_write_simple_molecule() {
        let mut mol = Molecule::new("methane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.6289, 0.6289, 0.6289));
        mol.atoms.push(Atom::new(2, "H", -0.6289, -0.6289, 0.6289));
        mol.atoms.push(Atom::new(3, "H", -0.6289, 0.6289, -0.6289));
        mol.atoms.push(Atom::new(4, "H", 0.6289, -0.6289, -0.6289));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let output = write_mol2_string(&mol).unwrap();

        assert!(output.contains("@<TRIPOS>MOLECULE"));
        assert!(output.contains("methane"));
        assert!(output.contains("@<TRIPOS>ATOM"));
        assert!(output.contains("@<TRIPOS>BOND"));
        assert!(output.contains("5 4 0 0 0")); // 5 atoms, 4 bonds
    }

    #[test]
    fn test_round_trip() {
        let mut mol = Molecule::new("test_molecule");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        let mol2_string = write_mol2_string(&mol).unwrap();
        let parsed = parse_mol2_string(&mol2_string).unwrap();

        assert_eq!(parsed.name, mol.name);
        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.bond_count(), mol.bond_count());

        // Check coordinates (with some tolerance for formatting)
        assert!((parsed.atoms[0].x - mol.atoms[0].x).abs() < 0.001);
        assert!((parsed.atoms[1].x - mol.atoms[1].x).abs() < 0.001);

        // Check bond order
        assert_eq!(parsed.bonds[0].order, mol.bonds[0].order);
    }

    #[test]
    fn test_write_aromatic_bonds() {
        let mut mol = Molecule::new("benzene");
        // Simple benzene ring (6 carbons)
        for i in 0..6 {
            let angle = std::f64::consts::PI * 2.0 * (i as f64) / 6.0;
            let x = 1.4 * angle.cos();
            let y = 1.4 * angle.sin();
            mol.atoms.push(Atom::new(i, "C", x, y, 0.0));
        }
        // Aromatic bonds
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }

        let output = write_mol2_string(&mol).unwrap();

        // All bonds should be written as "ar"
        assert!(output.contains(" ar"));
        // Atoms should have aromatic type
        assert!(output.contains("C.ar"));
    }

    #[test]
    fn test_bond_type_conversion() {
        assert_eq!(order_to_mol2_bond_type(BondOrder::Single), "1");
        assert_eq!(order_to_mol2_bond_type(BondOrder::Double), "2");
        assert_eq!(order_to_mol2_bond_type(BondOrder::Triple), "3");
        assert_eq!(order_to_mol2_bond_type(BondOrder::Aromatic), "ar");
    }

    #[test]
    fn test_write_multi() {
        let mut mol1 = Molecule::new("mol1");
        mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

        let mut mol2 = Molecule::new("mol2");
        mol2.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));

        let mut buffer = Vec::new();
        write_mol2_multi(&mut buffer, &[mol1, mol2]).unwrap();
        let output = String::from_utf8_lossy(&buffer);

        // Should contain two MOLECULE sections
        assert_eq!(output.matches("@<TRIPOS>MOLECULE").count(), 2);
        assert!(output.contains("mol1"));
        assert!(output.contains("mol2"));
    }
}
