use std::io::Write;

use crate::error::Result;
use crate::molecule::Molecule;

/// Writes a molecule to SDF V2000 format.
pub fn write_sdf<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
    // Line 1: Molecule name
    writeln!(writer, "{}", molecule.name)?;

    // Line 2: Program/timestamp line
    if let Some(ref prog_line) = molecule.program_line {
        writeln!(writer, "{}", prog_line)?;
    } else {
        writeln!(writer)?;
    }

    // Line 3: Comment
    if let Some(ref comment) = molecule.comment {
        writeln!(writer, "{}", comment)?;
    } else {
        writeln!(writer)?;
    }

    // Line 4: Counts line
    writeln!(
        writer,
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000",
        molecule.atom_count(),
        molecule.bond_count()
    )?;

    // Atom block
    for atom in &molecule.atoms {
        // Charge code mapping (reverse of parsing)
        let charge_code: u8 = match atom.formal_charge {
            3 => 1,
            2 => 2,
            1 => 3,
            0 => 0,
            -1 => 5,
            -2 => 6,
            -3 => 7,
            _ => 0, // Larger charges need M  CHG line
        };

        writeln!(
            writer,
            "{:10.4}{:10.4}{:10.4} {:3}{:2}{:3}  0  0  0  0  0  0  0  0  0  0",
            atom.x, atom.y, atom.z, atom.element, atom.mass_difference, charge_code
        )?;
    }

    // Bond block
    for bond in &molecule.bonds {
        writeln!(
            writer,
            "{:3}{:3}{:3}{:3}  0  0  0",
            bond.atom1 + 1, // Convert to 1-based
            bond.atom2 + 1,
            bond.order.to_sdf(),
            bond.stereo.to_sdf()
        )?;
    }

    // Write M  CHG lines for charges outside -3 to +3 range
    let charged_atoms: Vec<_> = molecule
        .atoms
        .iter()
        .filter(|a| a.formal_charge < -3 || a.formal_charge > 3)
        .collect();

    if !charged_atoms.is_empty() {
        for chunk in charged_atoms.chunks(8) {
            write!(writer, "M  CHG{:3}", chunk.len())?;
            for atom in chunk {
                write!(writer, " {:3}{:4}", atom.index + 1, atom.formal_charge)?;
            }
            writeln!(writer)?;
        }
    }

    // Write M  ISO lines for isotopes
    let isotope_atoms: Vec<_> = molecule
        .atoms
        .iter()
        .filter(|a| a.mass_difference != 0)
        .collect();

    if !isotope_atoms.is_empty() {
        for chunk in isotope_atoms.chunks(8) {
            write!(writer, "M  ISO{:3}", chunk.len())?;
            for atom in chunk {
                write!(writer, " {:3}{:4}", atom.index + 1, atom.mass_difference)?;
            }
            writeln!(writer)?;
        }
    }

    writeln!(writer, "M  END")?;

    // Data block (properties)
    for (key, value) in &molecule.properties {
        writeln!(writer, "> <{}>", key)?;
        writeln!(writer, "{}", value)?;
        writeln!(writer)?;
    }

    // Record separator
    writeln!(writer, "$$$$")?;

    Ok(())
}

/// Writes a molecule to an SDF string.
pub fn write_sdf_string(molecule: &Molecule) -> Result<String> {
    let mut buffer = Vec::new();
    write_sdf(&mut buffer, molecule)?;
    Ok(String::from_utf8_lossy(&buffer).to_string())
}

/// Writes multiple molecules to SDF format.
pub fn write_sdf_multi<W: Write>(writer: &mut W, molecules: &[Molecule]) -> Result<()> {
    for mol in molecules {
        write_sdf(writer, mol)?;
    }
    Ok(())
}

/// Writes a molecule to an SDF file.
pub fn write_sdf_file<P: AsRef<std::path::Path>>(path: P, molecule: &Molecule) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_sdf(&mut writer, molecule)
}

/// Writes multiple molecules to an SDF file.
pub fn write_sdf_file_multi<P: AsRef<std::path::Path>>(
    path: P,
    molecules: &[Molecule],
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_sdf_multi(&mut writer, molecules)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};
    use crate::parser::parse_sdf_string;

    #[test]
    fn test_round_trip() {
        let mut mol = Molecule::new("test_molecule");

        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        mol.set_property("TEST_PROP", "test_value");

        let sdf_string = write_sdf_string(&mol).unwrap();
        let parsed = parse_sdf_string(&sdf_string).unwrap();

        assert_eq!(parsed.name, mol.name);
        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.bond_count(), mol.bond_count());
        assert_eq!(parsed.get_property("TEST_PROP"), Some("test_value"));

        // Check coordinates
        assert!((parsed.atoms[0].x - mol.atoms[0].x).abs() < 0.001);
        assert!((parsed.atoms[1].x - mol.atoms[1].x).abs() < 0.001);

        // Check bond order
        assert_eq!(parsed.bonds[0].order, mol.bonds[0].order);
    }
}
