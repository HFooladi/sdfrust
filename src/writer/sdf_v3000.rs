//! SDF V3000 format writer.
//!
//! Writes molecules in V3000 format, which supports:
//! - Molecules with >999 atoms/bonds
//! - Extended bond types (coordination, hydrogen)
//! - Enhanced stereochemistry
//! - SGroups and collections

use std::io::Write;

use crate::error::Result;
use crate::molecule::{Molecule, SdfFormat};
use crate::stereogroup::StereoGroupType;

/// Writes a molecule to SDF V3000 format.
pub fn write_sdf_v3000<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
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

    // Line 4: Counts line (V3000 indicator)
    writeln!(writer, "  0  0  0     0  0            999 V3000")?;

    // Begin CTAB block
    writeln!(writer, "M  V30 BEGIN CTAB")?;

    // COUNTS line
    let sgroup_count = molecule.sgroups.len();
    writeln!(
        writer,
        "M  V30 COUNTS {} {} {} 0 0",
        molecule.atom_count(),
        molecule.bond_count(),
        sgroup_count
    )?;

    // ATOM block
    writeln!(writer, "M  V30 BEGIN ATOM")?;
    for (i, atom) in molecule.atoms.iter().enumerate() {
        let atom_id = atom.v3000_id.unwrap_or((i + 1) as u32);
        let aamap = atom.atom_atom_mapping.unwrap_or(0);

        // Build optional key=value pairs
        let mut extras = String::new();

        if atom.formal_charge != 0 {
            extras.push_str(&format!(" CHG={}", atom.formal_charge));
        }

        if let Some(rad) = atom.radical {
            if rad > 0 {
                extras.push_str(&format!(" RAD={}", rad));
            }
        }

        if atom.mass_difference != 0 {
            extras.push_str(&format!(" MASS={}", atom.mass_difference));
        }

        if let Some(val) = atom.valence {
            if val > 0 {
                extras.push_str(&format!(" VAL={}", val));
            }
        }

        if let Some(hcount) = atom.hydrogen_count {
            if hcount > 0 {
                extras.push_str(&format!(" HCOUNT={}", hcount));
            }
        }

        if let Some(cfg) = atom.stereo_parity {
            if cfg > 0 {
                extras.push_str(&format!(" CFG={}", cfg));
            }
        }

        if let Some(rg) = atom.rgroup_label {
            extras.push_str(&format!(" RGROUPS=(1 {})", rg));
        }

        writeln!(
            writer,
            "M  V30 {} {} {:.4} {:.4} {:.4} {}{}",
            atom_id, atom.element, atom.x, atom.y, atom.z, aamap, extras
        )?;
    }
    writeln!(writer, "M  V30 END ATOM")?;

    // BOND block
    writeln!(writer, "M  V30 BEGIN BOND")?;
    for (i, bond) in molecule.bonds.iter().enumerate() {
        let bond_id = bond.v3000_id.unwrap_or((i + 1) as u32);

        // Get atom IDs (use stored V3000 IDs or compute from indices)
        let atom1_id = molecule.atoms[bond.atom1]
            .v3000_id
            .unwrap_or((bond.atom1 + 1) as u32);
        let atom2_id = molecule.atoms[bond.atom2]
            .v3000_id
            .unwrap_or((bond.atom2 + 1) as u32);

        let bond_type = bond.order.to_sdf();

        // Build optional key=value pairs
        let mut extras = String::new();

        // Stereo configuration
        let cfg = match bond.stereo {
            crate::bond::BondStereo::Up => 1,
            crate::bond::BondStereo::Down => 3,
            crate::bond::BondStereo::Either => 2,
            crate::bond::BondStereo::None => 0,
        };
        if cfg > 0 {
            extras.push_str(&format!(" CFG={}", cfg));
        }

        if let Some(topo) = bond.topology {
            if topo > 0 {
                extras.push_str(&format!(" TOPO={}", topo));
            }
        }

        if let Some(rxctr) = bond.reacting_center {
            if rxctr > 0 {
                extras.push_str(&format!(" RXCTR={}", rxctr));
            }
        }

        writeln!(
            writer,
            "M  V30 {} {} {} {}{}",
            bond_id, bond_type, atom1_id, atom2_id, extras
        )?;
    }
    writeln!(writer, "M  V30 END BOND")?;

    // SGROUP block (if any)
    if !molecule.sgroups.is_empty() {
        writeln!(writer, "M  V30 BEGIN SGROUP")?;
        for sgroup in &molecule.sgroups {
            let type_str = sgroup.sgroup_type.to_v3000_str();

            let mut line = format!("M  V30 {} {}", sgroup.id, type_str);

            // Add atoms
            if !sgroup.atoms.is_empty() {
                line.push_str(&format!(" ATOMS=({}", sgroup.atoms.len()));
                for &idx in &sgroup.atoms {
                    let atom_id = molecule.atoms[idx].v3000_id.unwrap_or((idx + 1) as u32);
                    line.push_str(&format!(" {}", atom_id));
                }
                line.push(')');
            }

            // Add label
            if let Some(ref label) = sgroup.label {
                line.push_str(&format!(" LABEL=\"{}\"", label));
            }

            // Add subscript
            if let Some(ref sub) = sgroup.subscript {
                line.push_str(&format!(" SUBSCRIPT=\"{}\"", sub));
            }

            // Add connectivity
            if let Some(ref conn) = sgroup.connectivity {
                line.push_str(&format!(" CONNECT={}", conn));
            }

            writeln!(writer, "{}", line)?;
        }
        writeln!(writer, "M  V30 END SGROUP")?;
    }

    // Enhanced stereochemistry (COLLECTION block)
    if !molecule.stereogroups.is_empty() {
        writeln!(writer, "M  V30 BEGIN COLLECTION")?;
        for sg in &molecule.stereogroups {
            let type_str = match sg.group_type {
                StereoGroupType::Absolute => "MDLV30/STEABS".to_string(),
                StereoGroupType::Or => format!("MDLV30/STEREL{}", sg.group_number),
                StereoGroupType::And => format!("MDLV30/STERAC{}", sg.group_number),
            };

            let mut line = format!("M  V30 {}", type_str);

            if !sg.atoms.is_empty() {
                line.push_str(&format!(" ATOMS=({}", sg.atoms.len()));
                for &idx in &sg.atoms {
                    let atom_id = molecule.atoms[idx].v3000_id.unwrap_or((idx + 1) as u32);
                    line.push_str(&format!(" {}", atom_id));
                }
                line.push(')');
            }

            writeln!(writer, "{}", line)?;
        }
        writeln!(writer, "M  V30 END COLLECTION")?;
    }

    // End CTAB block
    writeln!(writer, "M  V30 END CTAB")?;

    // M  END
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

/// Writes a molecule to a V3000 SDF string.
pub fn write_sdf_v3000_string(molecule: &Molecule) -> Result<String> {
    let mut buffer = Vec::new();
    write_sdf_v3000(&mut buffer, molecule)?;
    Ok(String::from_utf8_lossy(&buffer).to_string())
}

/// Writes multiple molecules to V3000 SDF format.
pub fn write_sdf_v3000_multi<W: Write>(writer: &mut W, molecules: &[Molecule]) -> Result<()> {
    for mol in molecules {
        write_sdf_v3000(writer, mol)?;
    }
    Ok(())
}

/// Writes a molecule to a V3000 SDF file.
pub fn write_sdf_v3000_file<P: AsRef<std::path::Path>>(path: P, molecule: &Molecule) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_sdf_v3000(&mut writer, molecule)
}

/// Writes multiple molecules to a V3000 SDF file.
pub fn write_sdf_v3000_file_multi<P: AsRef<std::path::Path>>(
    path: P,
    molecules: &[Molecule],
) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_sdf_v3000_multi(&mut writer, molecules)
}

/// Writes a molecule to the appropriate SDF format (auto-detect).
///
/// Uses V3000 format if the molecule requires it (>999 atoms/bonds,
/// V3000-specific features), otherwise uses V2000.
pub fn write_sdf_auto<W: Write>(writer: &mut W, molecule: &Molecule) -> Result<()> {
    if molecule.needs_v3000() || molecule.format_version == SdfFormat::V3000 {
        write_sdf_v3000(writer, molecule)
    } else {
        crate::writer::write_sdf(writer, molecule)
    }
}

/// Writes a molecule to SDF string with automatic format selection.
pub fn write_sdf_auto_string(molecule: &Molecule) -> Result<String> {
    let mut buffer = Vec::new();
    write_sdf_auto(&mut buffer, molecule)?;
    Ok(String::from_utf8_lossy(&buffer).to_string())
}

/// Writes a molecule to SDF file with automatic format selection.
pub fn write_sdf_auto_file<P: AsRef<std::path::Path>>(path: P, molecule: &Molecule) -> Result<()> {
    let file = std::fs::File::create(path)?;
    let mut writer = std::io::BufWriter::new(file);
    write_sdf_auto(&mut writer, molecule)
}

/// Returns true if the molecule needs V3000 format.
///
/// A molecule needs V3000 format if:
/// - It has more than 999 atoms or bonds
/// - It has V3000-specific features (stereogroups, sgroups, collections)
/// - It has extended bond types (coordination, hydrogen)
/// - It was originally parsed from V3000 format
pub fn needs_v3000(molecule: &Molecule) -> bool {
    molecule.needs_v3000()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};
    use crate::parser::parse_sdf_v3000_string;

    #[test]
    fn test_write_v3000_simple() {
        let mut mol = Molecule::new("test_molecule");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        let sdf_string = write_sdf_v3000_string(&mol).unwrap();

        assert!(sdf_string.contains("V3000"));
        assert!(sdf_string.contains("COUNTS 2 1"));
        assert!(sdf_string.contains("BEGIN ATOM"));
        assert!(sdf_string.contains("END ATOM"));
        assert!(sdf_string.contains("BEGIN BOND"));
        assert!(sdf_string.contains("END BOND"));
    }

    #[test]
    fn test_v3000_round_trip() {
        let mut mol = Molecule::new("round_trip_test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "N", 1.5, 0.0, 0.0));
        mol.atoms[1].formal_charge = 1;
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.set_property("TEST_PROP", "test_value");

        let sdf_string = write_sdf_v3000_string(&mol).unwrap();
        let parsed = parse_sdf_v3000_string(&sdf_string).unwrap();

        assert_eq!(parsed.name, mol.name);
        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.bond_count(), mol.bond_count());
        assert_eq!(parsed.atoms[1].formal_charge, 1);
        assert_eq!(parsed.get_property("TEST_PROP"), Some("test_value"));
    }

    #[test]
    fn test_write_v3000_with_charge() {
        let mut mol = Molecule::new("charged");
        mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
        mol.atoms[0].formal_charge = 1;
        mol.atoms.push(Atom::new(1, "O", 1.5, 0.0, 0.0));
        mol.atoms[1].formal_charge = -1;
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let sdf_string = write_sdf_v3000_string(&mol).unwrap();

        assert!(sdf_string.contains("CHG=1"));
        assert!(sdf_string.contains("CHG=-1"));
    }

    #[test]
    fn test_needs_v3000() {
        let mut mol = Molecule::new("small");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 0, BondOrder::Single));

        assert!(!needs_v3000(&mol));

        // Add coordination bond
        mol.bonds[0].order = BondOrder::Coordination;
        assert!(needs_v3000(&mol));
    }

    #[test]
    fn test_write_sdf_auto() {
        // Small molecule should use V2000
        let mut mol = Molecule::new("small");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 0, BondOrder::Single));

        let sdf_string = write_sdf_auto_string(&mol).unwrap();
        assert!(sdf_string.contains("V2000"));
        assert!(!sdf_string.contains("V3000"));

        // Molecule with coordination bond should use V3000
        mol.bonds[0].order = BondOrder::Coordination;
        let sdf_string = write_sdf_auto_string(&mol).unwrap();
        assert!(sdf_string.contains("V3000"));
    }
}
