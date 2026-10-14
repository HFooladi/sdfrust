//! Regression tests for correctness bugs found in the v0.6.0 audit.
//!
//! Each test pins down a case where sdfrust previously returned wrong data,
//! panicked, or disagreed with RDKit/OpenBabel on the same input.

use std::io::Cursor;

use sdfrust::featurize::ogb;
use sdfrust::{
    Atom, Bond, BondOrder, Molecule, SdfError, SdfIterator, parse_auto_string_multi,
    parse_mol2_string, parse_sdf_string, parse_sdf_string_multi, parse_sdf_v3000_string,
    parse_xyz_string, write_sdf_auto_string, write_sdf_string, write_sdf_v3000_string,
};

const V3000_BENZENE: &str = include_str!("test_data/v3000_benzene.sdf");
const ASPIRIN: &str = include_str!("test_data/aspirin.sdf");

fn methane_block(name: &str) -> String {
    format!(
        "{name}\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n"
    )
}

fn carbon_chain(n: usize) -> Molecule {
    let mut mol = Molecule::new("chain");
    for i in 0..n {
        mol.atoms.push(Atom::new(i, "C", i as f64 * 1.5, 0.0, 0.0));
    }
    for i in 1..n {
        mol.bonds.push(Bond::new(i - 1, i, BondOrder::Single));
    }
    mol
}

// ---------------------------------------------------------------------------
// Record framing
// ---------------------------------------------------------------------------

#[test]
fn v3000_record_through_v2000_entry_point_is_parsed() {
    let mol = parse_sdf_string(V3000_BENZENE).unwrap();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
    assert_eq!(mol.get_property("FORMULA"), Some("C6H6"));
}

#[test]
fn mixed_v2000_v3000_file_parses_every_record() {
    let content = format!("{ASPIRIN}{V3000_BENZENE}{}", methane_block("methane"));
    let mols = parse_sdf_string_multi(&content).unwrap();
    assert_eq!(
        mols.iter().map(|m| m.atom_count()).collect::<Vec<_>>(),
        vec![21, 6, 1]
    );
    let auto = parse_auto_string_multi(&content).unwrap();
    assert_eq!(auto.len(), 3);
}

#[test]
fn trailing_blank_lines_after_last_record_are_ignored() {
    let content = format!("{}\n\n  \n", methane_block("a"));
    assert_eq!(parse_sdf_string_multi(&content).unwrap().len(), 1);

    let v3000 = format!("{V3000_BENZENE}\n\n");
    assert_eq!(
        sdfrust::parse_sdf_v3000_string_multi(&v3000).unwrap().len(),
        1
    );
}

#[test]
fn record_with_blank_header_lines_still_parses() {
    let mols = parse_sdf_string_multi(&methane_block("")).unwrap();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].atom_count(), 1);
}

#[test]
fn iterator_resumes_after_malformed_record() {
    let bad = "bad\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\nGARBAGE\nM  END\n$$$$\n";
    let content = format!("{}{bad}{}", methane_block("first"), methane_block("third"));
    let results: Vec<_> = SdfIterator::new(Cursor::new(content)).collect();
    assert_eq!(results.len(), 3);
    assert!(results[0].is_ok());
    assert!(results[1].is_err());
    assert_eq!(results[2].as_ref().unwrap().name, "third");
}

#[test]
fn missing_m_end_does_not_swallow_next_record() {
    let truncated = "first\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0\n$$$$\n";
    let content = format!("{truncated}{}", methane_block("second"));
    let mols = parse_sdf_string_multi(&content).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[1].name, "second");
}

#[test]
fn non_ascii_in_fixed_width_line_is_an_error_not_a_panic() {
    let content = "x\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0µÅ00 C   0  0\nM  END\n$$$$\n";
    assert!(parse_sdf_string(content).is_err());
}

#[test]
fn invalid_utf8_in_data_block_is_replaced_not_fatal() {
    let mut bytes = methane_block("m").into_bytes();
    let tail = bytes.split_off(bytes.len() - "$$$$\n".len());
    bytes.extend_from_slice(b"> <NOTE>\ncaf\xe9\n\n");
    bytes.extend_from_slice(&tail);
    let mol = sdfrust::SdfParser::new(Cursor::new(bytes))
        .parse_molecule()
        .unwrap()
        .unwrap();
    assert_eq!(mol.get_property("NOTE"), Some("caf\u{FFFD}"));
}

#[test]
fn data_item_ends_at_blank_line() {
    let mut content = methane_block("m");
    content.truncate(content.len() - "$$$$\n".len());
    content.push_str("> <A>\nline1\nline2\n\nstray text\n> <B>\nb\n\n$$$$\n");
    let mol = parse_sdf_string(&content).unwrap();
    assert_eq!(mol.get_property("A"), Some("line1\nline2"));
    assert_eq!(mol.get_property("B"), Some("b"));
}

#[test]
fn huge_declared_counts_are_errors_not_aborts() {
    let v3000 = "big\n\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 999999999999 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n";
    assert!(parse_sdf_v3000_string(v3000).is_err());
    assert!(parse_xyz_string("999999999999\nbig\nC 0 0 0\n").is_err());
}

// ---------------------------------------------------------------------------
// Charges, isotopes, radicals
// ---------------------------------------------------------------------------

#[test]
fn m_chg_supersedes_atom_block_charges() {
    // Atom block says +3 (code 1) on atom 1; M  CHG only mentions atom 2.
    let content = "x\n\n\n  2  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 N   0  1  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0\nM  CHG  1   2  -1\nM  END\n$$$$\n";
    let mol = parse_sdf_string(content).unwrap();
    assert_eq!(mol.atoms[0].formal_charge, 0);
    assert_eq!(mol.atoms[1].formal_charge, -1);
}

#[test]
fn small_and_large_charges_survive_round_trip() {
    let mut mol = carbon_chain(3);
    mol.atoms[0].formal_charge = 1;
    mol.atoms[1].formal_charge = 4;
    mol.atoms[2].formal_charge = -2;
    let text = write_sdf_string(&mol).unwrap();
    assert!(text.contains("M  CHG  3   1   1   2   4   3  -2"));
    let back = parse_sdf_string(&text).unwrap();
    let charges: Vec<i8> = back.atoms.iter().map(|a| a.formal_charge).collect();
    assert_eq!(charges, vec![1, 4, -2]);
}

#[test]
fn isotopes_are_absolute_mass_numbers_and_round_trip() {
    for mass in [2u16, 13, 235] {
        let content = format!(
            "x\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  ISO  1   1 {mass:3}\nM  END\n$$$$\n"
        );
        let mol = parse_sdf_string(&content).unwrap();
        assert_eq!(mol.atoms[0].isotope, Some(mass));
        assert_eq!(mol.atoms[0].mass_difference, 0);

        let back = parse_sdf_string(&write_sdf_string(&mol).unwrap()).unwrap();
        assert_eq!(back.atoms[0].isotope, Some(mass));

        let v3 = parse_sdf_v3000_string(&write_sdf_v3000_string(&mol).unwrap()).unwrap();
        assert_eq!(v3.atoms[0].isotope, Some(mass));
    }
}

#[test]
fn v3000_decimal_mass_is_an_isotope() {
    let content = V3000_BENZENE.replace(
        "M  V30 1 C 0.0000 1.4000 0.0000 0",
        "M  V30 1 C 0.0000 1.4000 0.0000 0 MASS=13.003",
    );
    let mol = parse_sdf_string(&content).unwrap();
    assert_eq!(mol.atoms[0].isotope, Some(13));
}

#[test]
fn radicals_are_parsed_and_written() {
    // Charge code 4 in the atom block is a doublet radical
    let code4 = "x\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  4  0  0  0  0\nM  END\n$$$$\n";
    let mol = parse_sdf_string(code4).unwrap();
    assert_eq!(mol.atoms[0].radical, Some(2));
    assert_eq!(mol.atoms[0].radical_electrons(), 1);
    assert_eq!(mol.atoms[0].formal_charge, 0);

    let text = write_sdf_string(&mol).unwrap();
    assert!(text.contains("M  RAD  1   1   2"));
    assert_eq!(parse_sdf_string(&text).unwrap().atoms[0].radical, Some(2));
    // Radicals no longer force V3000 output
    assert!(write_sdf_auto_string(&mol).unwrap().contains("V2000"));
}

// ---------------------------------------------------------------------------
// Writer limits and validation
// ---------------------------------------------------------------------------

#[test]
fn v2000_writer_rejects_more_than_999_atoms() {
    let mol = carbon_chain(1000);
    assert!(matches!(
        write_sdf_string(&mol),
        Err(SdfError::Unwritable(_))
    ));
    let auto = write_sdf_auto_string(&mol).unwrap();
    assert!(auto.contains("V3000"));
    assert_eq!(parse_sdf_string(&auto).unwrap().atom_count(), 1000);
}

#[test]
fn v2000_writer_rejects_coordinates_too_wide_for_columns() {
    let mut mol = carbon_chain(1);
    mol.atoms[0].x = -12345.0;
    assert!(matches!(
        write_sdf_string(&mol),
        Err(SdfError::Unwritable(_))
    ));
    assert!(write_sdf_auto_string(&mol).unwrap().contains("V3000"));
}

#[test]
fn writers_reject_bonds_to_missing_atoms() {
    let mut mol = carbon_chain(1);
    mol.bonds.push(Bond::new(0, 7, BondOrder::Single));
    assert!(matches!(
        write_sdf_string(&mol),
        Err(SdfError::InvalidAtomIndex { .. })
    ));
    assert!(matches!(
        write_sdf_v3000_string(&mol),
        Err(SdfError::InvalidAtomIndex { .. })
    ));
}

#[test]
fn writer_uses_atom_position_not_stale_index_field() {
    let mut mol = carbon_chain(2);
    mol.atoms[1].index = 42; // stale/user-assigned index
    mol.atoms[1].formal_charge = -1;
    let back = parse_sdf_string(&write_sdf_string(&mol).unwrap()).unwrap();
    assert_eq!(back.atoms[1].formal_charge, -1);
}

#[test]
fn property_output_is_sorted_and_deterministic() {
    let mut mol = carbon_chain(1);
    for key in ["zeta", "alpha", "mid"] {
        mol.set_property(key, key);
    }
    let text = write_sdf_string(&mol).unwrap();
    let a = text.find("<alpha>").unwrap();
    let m = text.find("<mid>").unwrap();
    let z = text.find("<zeta>").unwrap();
    assert!(a < m && m < z);
}

// ---------------------------------------------------------------------------
// MOL2 / XYZ
// ---------------------------------------------------------------------------

const MOL2_AMIDE: &str = "@<TRIPOS>MOLECULE
acetamide
 5 4 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
     10 C1  0.0 0.0 0.0 C.3  1 LIG  -0.20
     20 C2  1.5 0.0 0.0 C.2  1 LIG   0.62
     30 O1  2.1 1.0 0.0 O.2  1 LIG  -0.55
     40 N1  2.1 -1.2 0.0 N.am 1 LIG  -0.70
     50 CL1 -1.0 0.0 0.0 CL  1 LIG  -0.10
@<TRIPOS>BOND
     1 10 20 1
     2 20 30 2
     3 20 40 am
     4 10 50 1
";

#[test]
fn mol2_amide_bond_is_single() {
    let mol = parse_mol2_string(MOL2_AMIDE).unwrap();
    assert_eq!(mol.bonds[2].order, BondOrder::Single);
}

#[test]
fn mol2_non_contiguous_atom_ids_are_resolved() {
    let mol = parse_mol2_string(MOL2_AMIDE).unwrap();
    let pairs: Vec<_> = mol.bonds.iter().map(|b| (b.atom1, b.atom2)).collect();
    assert_eq!(pairs, vec![(0, 1), (1, 2), (1, 3), (0, 4)]);
}

#[test]
fn mol2_partial_charges_and_element_case() {
    let mol = parse_mol2_string(MOL2_AMIDE).unwrap();
    assert_eq!(mol.atoms[1].partial_charge, Some(0.62));
    assert!(mol.atoms.iter().all(|a| a.formal_charge == 0));
    assert_eq!(mol.atoms[4].element, "Cl");

    let no_charges = MOL2_AMIDE.replace("USER_CHARGES", "NO_CHARGES");
    let mol = parse_mol2_string(&no_charges).unwrap();
    assert!(mol.atoms.iter().all(|a| a.partial_charge.is_none()));
}

#[test]
fn mol2_quaternary_nitrogen_is_positive() {
    let mol2 = "@<TRIPOS>MOLECULE\nnh4\n 1 0 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n 1 N1 0 0 0 N.4 1 LIG\n";
    assert_eq!(parse_mol2_string(mol2).unwrap().atoms[0].formal_charge, 1);
}

#[test]
fn xyz_supports_all_atomic_numbers() {
    let mol = parse_xyz_string("3\nmetals\n78 0 0 0\n46 1 0 0\n118 2 0 0\n").unwrap();
    let elements: Vec<_> = mol.atoms.iter().map(|a| a.element.as_str()).collect();
    assert_eq!(elements, vec!["Pt", "Pd", "Og"]);
}

// ---------------------------------------------------------------------------
// OGB encoding (indices into ogb.utils.features.allowable_features)
// ---------------------------------------------------------------------------

#[test]
fn ogb_atom_features_match_ogb_reference_for_aspirin() {
    // Reference vectors computed with ogb 1.3.6 atom_to_feature_vector on RDKit's
    // parse of aspirin.sdf (removeHs=False).
    let mol = parse_sdf_string(ASPIRIN).unwrap();
    let feats = ogb::ogb_atom_features(&mol).features;
    // O (hydroxyl), C (aromatic), H
    let o_idx = mol.atoms.iter().position(|a| a.element == "O").unwrap();
    assert_eq!(feats[o_idx][0], 7);
    let h_idx = mol.atoms.iter().position(|a| a.element == "H").unwrap();
    assert_eq!(feats[h_idx][0], 0);
    assert_eq!(feats[h_idx][6], 5); // S hybridization → misc
    let dims = [119, 5, 12, 12, 10, 6, 6, 2, 2];
    for f in &feats {
        for (v, d) in f.iter().zip(dims) {
            assert!((0..d).contains(v), "{f:?} out of OGB range");
        }
    }
}

// ---------------------------------------------------------------------------
// Phase 2: implicit hydrogens and masses (values from RDKit 2025.09)
// ---------------------------------------------------------------------------

fn heavy(element: &str, charge: i8) -> Molecule {
    let mut mol = Molecule::new(element);
    let mut atom = Atom::new(0, element, 0.0, 0.0, 0.0);
    atom.formal_charge = charge;
    mol.atoms.push(atom);
    mol
}

#[test]
fn implicit_hydrogens_follow_charge_aware_valence() {
    use sdfrust::descriptors::valence::implicit_hydrogen_count;
    for (element, charge, expected) in [
        ("C", 0, 4),
        ("N", 1, 4),  // ammonium
        ("O", -1, 1), // hydroxide
        ("O", 1, 3),  // hydronium
        ("C", 1, 3),  // methyl cation
        ("C", -1, 3), // methyl anion
        ("B", -1, 4), // borohydride
        ("Fe", 2, 0), // metals take no implicit H
        ("Cl", -1, 0),
    ] {
        assert_eq!(
            implicit_hydrogen_count(&heavy(element, charge), 0),
            expected,
            "{element}{charge:+}"
        );
    }
}

#[test]
fn methyl_radical_from_molfile_keeps_three_hydrogens() {
    // As written by RDKit for [CH3]: valence field 3 and M  RAD doublet
    let block = "\n     RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  3  0  0  0  0  0  0\nM  RAD  1   1   2\nM  END\n";
    let mol = parse_sdf_string(block).unwrap();
    assert_eq!(mol.formula(), "CH3");
}

#[test]
fn masses_and_formula_include_implicit_h_isotopes_and_charge() {
    // Acetate, heavy atoms only: C2H3O2-
    let mut mol = carbon_chain(2);
    mol.atoms.push(Atom::new(2, "O", 0.0, 1.0, 0.0));
    mol.atoms.push(Atom::new(3, "O", 0.0, -1.0, 0.0));
    mol.atoms[3].formal_charge = -1;
    mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
    mol.bonds.push(Bond::new(1, 3, BondOrder::Single));
    assert_eq!(mol.formula(), "C2H3O2-");
    assert!((mol.molecular_weight().unwrap() - 59.044).abs() < 1e-3);
    assert!((mol.exact_mass().unwrap() - 59.013853).abs() < 1e-5);

    // 13CH4 uses the isotope mass
    let mut c13 = heavy("C", 0);
    c13.atoms[0].isotope = Some(13);
    assert!((c13.molecular_weight().unwrap() - 17.0354).abs() < 1e-3);
    assert_eq!(c13.formula(), "CH4");

    // Uranium is now a known element
    assert!(heavy("U", 0).molecular_weight().is_some());
}
