use std::io::Write;

use crate::bond::BondStereo;
use crate::error::{Result, SdfError};
use crate::molecule::Molecule;

/// Largest absolute coordinate that fits the V2000 `%10.4f` fields.
const V2000_MAX_COORD: f64 = 9999.9999;

/// Checks that every bond (and stereo/sgroup reference) points at an existing atom.
pub(crate) fn validate_atom_references(molecule: &Molecule) -> Result<()> {
    let n = molecule.atom_count();
    let check = |idx: usize| {
        if idx < n {
            Ok(())
        } else {
            Err(SdfError::InvalidAtomIndex {
                index: idx,
                atom_count: n,
            })
        }
    };
    for bond in &molecule.bonds {
        check(bond.atom1)?;
        check(bond.atom2)?;
    }
    for &idx in molecule
        .sgroups
        .iter()
        .flat_map(|g| g.atoms.iter())
        .chain(molecule.stereogroups.iter().flat_map(|g| g.atoms.iter()))
    {
        check(idx)?;
    }
    Ok(())
}

/// Returns why `molecule` cannot be represented in V2000, if it can't.
pub(crate) fn v2000_limitation(molecule: &Molecule) -> Option<String> {
    if molecule.atom_count() > 999 || molecule.bond_count() > 999 {
        return Some(format!(
            "{} atoms / {} bonds exceed the V2000 limit of 999",
            molecule.atom_count(),
            molecule.bond_count()
        ));
    }
    for (i, atom) in molecule.atoms.iter().enumerate() {
        if [atom.x, atom.y, atom.z]
            .iter()
            .any(|c| !c.is_finite() || c.abs() > V2000_MAX_COORD)
        {
            return Some(format!(
                "atom {} coordinates do not fit V2000 columns",
                i + 1
            ));
        }
        if atom.element.len() > 3 || !atom.element.is_ascii() {
            return Some(format!(
                "atom {} symbol '{}' does not fit V2000",
                i + 1,
                atom.element
            ));
        }
    }
    None
}

/// Writes a molecule to SDF V2000 format.
///
/// Returns [`SdfError::Unwritable`] if the molecule exceeds V2000 limits (more
/// than 999 atoms or bonds, coordinates beyond ±9999.9999, long atom symbols);
/// use [`crate::write_sdf_auto`] to fall back to V3000 automatically.
///
/// Formal charges, isotopes and radicals are written to `M  CHG`, `M  ISO` and
/// `M  RAD` lines (which per the CTfile spec supersede the atom block); charges
/// of ±1..±3 are also written to the atom block for legacy readers.
pub fn write_sdf<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
    if let Some(reason) = v2000_limitation(molecule) {
        return Err(SdfError::Unwritable(reason));
    }
    validate_atom_references(molecule)?;

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
        // Charge code mapping (reverse of parsing); M  CHG below is authoritative
        let charge_code: u8 = match atom.formal_charge {
            3 => 1,
            2 => 2,
            1 => 3,
            -1 => 5,
            -2 => 6,
            -3 => 7,
            _ => 0,
        };
        // The atom-block mass difference is only meaningful without an explicit isotope
        let mass_diff = if atom.isotope.is_none() && (-3..=4).contains(&atom.mass_difference) {
            atom.mass_difference
        } else {
            0
        };

        writeln!(
            writer,
            "{:10.4}{:10.4}{:10.4} {:3}{:2}{:3}  0  0  0  0  0  0  0  0  0  0",
            atom.x, atom.y, atom.z, atom.element, mass_diff, charge_code
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
            if bond.order == crate::bond::BondOrder::Double && bond.stereo == BondStereo::Either {
                3 // V2000 "cis or trans" for double bonds
            } else {
                bond.stereo.to_sdf()
            }
        )?;
    }

    // Property lines (atom numbers are positions in the atom list, 1-based)
    let charges: Vec<(usize, i32)> = molecule
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.formal_charge != 0)
        .map(|(i, a)| (i + 1, a.formal_charge as i32))
        .collect();
    let radicals: Vec<(usize, i32)> = molecule
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| match a.radical {
            Some(r @ 1..=3) => Some((i + 1, r as i32)),
            _ => None,
        })
        .collect();
    let isotopes: Vec<(usize, i32)> = molecule
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| a.isotope.map(|m| (i + 1, m as i32)))
        .collect();

    write_property_lines(writer, "CHG", &charges)?;
    write_property_lines(writer, "RAD", &radicals)?;
    write_property_lines(writer, "ISO", &isotopes)?;

    writeln!(writer, "M  END")?;

    write_data_items(writer, molecule)?;

    // Record separator
    writeln!(writer, "$$$$")?;

    Ok(())
}

/// Writes `M  XXX` property lines, at most 8 entries per line.
fn write_property_lines<W: Write>(
    writer: &mut W,
    tag: &str,
    entries: &[(usize, i32)],
) -> Result<()> {
    for chunk in entries.chunks(8) {
        write!(writer, "M  {}{:3}", tag, chunk.len())?;
        for (atom, value) in chunk {
            write!(writer, " {:3} {:3}", atom, value)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Writes the `> <name>` data items, sorted by name so output is deterministic.
pub(crate) fn write_data_items<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
    let mut items: Vec<_> = molecule.properties.iter().collect();
    items.sort_by(|a, b| a.0.cmp(b.0));
    for (key, value) in items {
        writeln!(writer, "> <{}>", key)?;
        writeln!(writer, "{}", value)?;
        writeln!(writer)?;
    }
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
