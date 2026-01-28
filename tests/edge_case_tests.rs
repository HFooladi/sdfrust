//! Edge case tests for sdfrust
//!
//! Tests for unusual, malformed, or boundary condition inputs.

use sdfrust::{
    parse_sdf_string, write_sdf_string, Atom, Bond, BondOrder, BondStereo, Molecule,
};

// ============================================================================
// Empty and Minimal Molecules
// ============================================================================

#[test]
fn test_empty_molecule_no_atoms() {
    let sdf = r#"empty_molecule


  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "empty_molecule");
    assert_eq!(mol.atom_count(), 0);
    assert_eq!(mol.bond_count(), 0);
    assert!(mol.is_empty());
}

#[test]
fn test_single_atom_no_bonds() {
    let sdf = r#"single_atom


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atom_count(), 1);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(mol.formula(), "C");
}

#[test]
fn test_molecule_with_only_hydrogens() {
    let sdf = r#"H2


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7400    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.formula(), "H2");
    assert_eq!(mol.atom_count(), 2);
}

// ============================================================================
// Special Characters and Names
// ============================================================================

#[test]
fn test_molecule_name_with_spaces() {
    let sdf = r#"molecule with spaces in name


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "molecule with spaces in name");
}

#[test]
fn test_molecule_name_with_special_chars() {
    let sdf = r#"mol-123_test(v2)[alpha]


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "mol-123_test(v2)[alpha]");
}

#[test]
fn test_empty_molecule_name() {
    let sdf = r#"


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "");
}

#[test]
fn test_property_with_special_characters() {
    let sdf = r#"test


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <SMILES>
C=C(C)C(=O)O[C@H]1C[C@@H]2CC[C@H]1C2

> <NAME_WITH_UNICODE>
Aspirin (acetylsalicylic acid)

$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert!(mol.get_property("SMILES").is_some());
    assert!(mol.get_property("NAME_WITH_UNICODE").is_some());
}

// ============================================================================
// Coordinate Edge Cases
// ============================================================================

#[test]
fn test_large_coordinates() {
    let sdf = r#"large_coords


  2  1  0  0  0  0  0  0  0  0999 V2000
  999.9999  999.9999  999.9999 C   0  0  0  0  0  0  0  0  0  0  0  0
 -999.9999 -999.9999 -999.9999 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert!((mol.atoms[0].x - 999.9999).abs() < 0.001);
    assert!((mol.atoms[1].x - (-999.9999)).abs() < 0.001);
}

#[test]
fn test_zero_coordinates() {
    let sdf = r#"zero_coords


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].x, 0.0);
    assert_eq!(mol.atoms[0].y, 0.0);
    assert_eq!(mol.atoms[0].z, 0.0);
}

#[test]
fn test_negative_coordinates() {
    let sdf = r#"negative_coords


  1  0  0  0  0  0  0  0  0  0999 V2000
   -1.5000   -2.5000   -3.5000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert!((mol.atoms[0].x - (-1.5)).abs() < 0.001);
    assert!((mol.atoms[0].y - (-2.5)).abs() < 0.001);
    assert!((mol.atoms[0].z - (-3.5)).abs() < 0.001);
}

// ============================================================================
// Charge Edge Cases
// ============================================================================

#[test]
fn test_positive_charges() {
    let sdf = r#"cation


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].formal_charge, 1);
    assert_eq!(mol.total_charge(), 1);
}

#[test]
fn test_negative_charges() {
    let sdf = r#"anion


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].formal_charge, -1);
    assert_eq!(mol.total_charge(), -1);
}

#[test]
fn test_m_chg_line_charges() {
    let sdf = r#"zwitterion


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  2   1   1   2  -1
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].formal_charge, 1);
    assert_eq!(mol.atoms[1].formal_charge, -1);
    assert_eq!(mol.total_charge(), 0);
}

#[test]
fn test_large_charges_via_m_chg() {
    let sdf = r#"highly_charged


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   3
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].formal_charge, 3);
}

// ============================================================================
// Bond Edge Cases
// ============================================================================

#[test]
fn test_all_bond_orders() {
    let sdf = r#"all_bonds


  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  3  4  3  0  0  0  0
  4  5  4  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.bonds[0].order, BondOrder::Single);
    assert_eq!(mol.bonds[1].order, BondOrder::Double);
    assert_eq!(mol.bonds[2].order, BondOrder::Triple);
    assert_eq!(mol.bonds[3].order, BondOrder::Aromatic);
}

#[test]
fn test_bond_stereo_types() {
    let sdf = r#"stereo_bonds


  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  2  3  1  6  0  0  0
  3  4  1  4  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.bonds[0].stereo, BondStereo::Up);
    assert_eq!(mol.bonds[1].stereo, BondStereo::Down);
    assert_eq!(mol.bonds[2].stereo, BondStereo::Either);
}

// ============================================================================
// Element Edge Cases
// ============================================================================

#[test]
fn test_two_letter_elements() {
    let sdf = r#"metals


  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Cu  0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Zn  0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 Mg  0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].element, "Fe");
    assert_eq!(mol.atoms[1].element, "Cu");
    assert_eq!(mol.atoms[2].element, "Zn");
    assert_eq!(mol.atoms[3].element, "Mg");
}

#[test]
fn test_halogens() {
    let sdf = r#"halogens


  4  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.formula(), "BrClFI");
}

#[test]
fn test_rare_elements() {
    let sdf = r#"rare


  3  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Se  0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Te  0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 As  0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.atoms[0].element, "Se");
    assert_eq!(mol.atoms[1].element, "Te");
    assert_eq!(mol.atoms[2].element, "As");
}

// ============================================================================
// Multi-molecule Edge Cases
// ============================================================================

#[test]
fn test_two_molecules_in_one_file() {
    let sdf = r#"mol1


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
mol2


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mols = sdfrust::parse_sdf_string_multi(sdf).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "mol1");
    assert_eq!(mols[1].name, "mol2");
    assert_eq!(mols[0].atoms[0].element, "C");
    assert_eq!(mols[1].atoms[0].element, "N");
}

#[test]
fn test_many_molecules() {
    let mut sdf = String::new();
    for i in 0..50 {
        sdf.push_str(&format!(
            r#"mol{}


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#,
            i
        ));
    }

    let mols = sdfrust::parse_sdf_string_multi(&sdf).unwrap();
    assert_eq!(mols.len(), 50);
    assert_eq!(mols[49].name, "mol49");
}

// ============================================================================
// Property Edge Cases
// ============================================================================

#[test]
fn test_empty_property_value() {
    let sdf = r#"test


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <EMPTY_PROP>

> <HAS_VALUE>
value

$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.get_property("EMPTY_PROP"), Some(""));
    assert_eq!(mol.get_property("HAS_VALUE"), Some("value"));
}

#[test]
fn test_multiline_property_value() {
    let sdf = r#"test


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <MULTILINE>
line1
line2
line3

$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    let value = mol.get_property("MULTILINE").unwrap();
    assert!(value.contains("line1"));
    assert!(value.contains("line2"));
    assert!(value.contains("line3"));
}

#[test]
fn test_many_properties() {
    let mut sdf = String::from(
        r#"test


  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"#,
    );

    for i in 0..100 {
        sdf.push_str(&format!(
            "> <PROP{}>\nvalue{}\n\n",
            i, i
        ));
    }
    sdf.push_str("$$$$\n");

    let mol = parse_sdf_string(&sdf).unwrap();
    assert_eq!(mol.properties.len(), 100);
    assert_eq!(mol.get_property("PROP0"), Some("value0"));
    assert_eq!(mol.get_property("PROP99"), Some("value99"));
}

// ============================================================================
// Round-trip Edge Cases
// ============================================================================

#[test]
fn test_roundtrip_empty_molecule() {
    let mol = Molecule::new("empty");
    let sdf = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&sdf).unwrap();
    assert_eq!(parsed.name, mol.name);
    assert_eq!(parsed.atom_count(), 0);
}

#[test]
fn test_roundtrip_with_charges() {
    let mut mol = Molecule::new("charged");
    let mut atom = Atom::new(0, "N", 0.0, 0.0, 0.0);
    atom.formal_charge = 1;
    mol.atoms.push(atom);

    let sdf = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&sdf).unwrap();
    assert_eq!(parsed.atoms[0].formal_charge, 1);
}

#[test]
fn test_roundtrip_with_stereo_bonds() {
    let mut mol = Molecule::new("stereo");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.bonds.push(Bond::with_stereo(0, 1, BondOrder::Single, BondStereo::Up));

    let sdf = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&sdf).unwrap();
    assert_eq!(parsed.bonds[0].stereo, BondStereo::Up);
}

#[test]
fn test_roundtrip_preserves_coordinates() {
    let mut mol = Molecule::new("coords");
    mol.atoms.push(Atom::new(0, "C", 1.2345, -6.7890, 0.1234));

    let sdf = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&sdf).unwrap();

    assert!((parsed.atoms[0].x - 1.2345).abs() < 0.001);
    assert!((parsed.atoms[0].y - (-6.7890)).abs() < 0.001);
    assert!((parsed.atoms[0].z - 0.1234).abs() < 0.001);
}

#[test]
fn test_roundtrip_preserves_properties() {
    let mut mol = Molecule::new("props");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.set_property("KEY1", "value1");
    mol.set_property("KEY2", "value2");

    let sdf = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&sdf).unwrap();

    assert_eq!(parsed.get_property("KEY1"), Some("value1"));
    assert_eq!(parsed.get_property("KEY2"), Some("value2"));
}

// ============================================================================
// Formula Edge Cases
// ============================================================================

#[test]
fn test_formula_single_element() {
    let mut mol = Molecule::new("single");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    assert_eq!(mol.formula(), "C");
}

#[test]
fn test_formula_ordering() {
    // Formula should be C, H, then alphabetical
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "N", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(3, "C", 0.0, 0.0, 0.0));

    assert_eq!(mol.formula(), "CHNO");
}

#[test]
fn test_formula_with_counts() {
    let mut mol = Molecule::new("test");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 6..12 {
        mol.atoms.push(Atom::new(i, "H", 0.0, 0.0, 0.0));
    }

    assert_eq!(mol.formula(), "C6H6");
}

// ============================================================================
// Geometry Edge Cases
// ============================================================================

#[test]
fn test_centroid_single_atom() {
    let mut mol = Molecule::new("single");
    mol.atoms.push(Atom::new(0, "C", 5.0, 3.0, 1.0));

    let (cx, cy, cz) = mol.centroid().unwrap();
    assert_eq!(cx, 5.0);
    assert_eq!(cy, 3.0);
    assert_eq!(cz, 1.0);
}

#[test]
fn test_centroid_empty_molecule() {
    let mol = Molecule::new("empty");
    assert!(mol.centroid().is_none());
}

#[test]
fn test_translate_empty_molecule() {
    let mut mol = Molecule::new("empty");
    mol.translate(1.0, 2.0, 3.0); // Should not panic
    assert!(mol.is_empty());
}

#[test]
fn test_center_empty_molecule() {
    let mut mol = Molecule::new("empty");
    mol.center(); // Should not panic
    assert!(mol.is_empty());
}

#[test]
fn test_atom_distance_same_position() {
    let a1 = Atom::new(0, "C", 1.0, 2.0, 3.0);
    let a2 = Atom::new(1, "C", 1.0, 2.0, 3.0);
    assert_eq!(a1.distance_to(&a2), 0.0);
}

// ============================================================================
// Error Handling Tests
// ============================================================================

#[test]
fn test_error_empty_file() {
    let result = parse_sdf_string("");
    assert!(result.is_err());
}

#[test]
fn test_error_invalid_atom_count() {
    let sdf = r#"test


  2  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;
    // Declares 2 atoms but only has 1
    let result = parse_sdf_string(sdf);
    assert!(result.is_err());
}

#[test]
fn test_error_invalid_bond_atom_index() {
    let sdf = r#"test


  1  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;
    // Bond references atom 5, but only 1 atom exists
    let result = parse_sdf_string(sdf);
    assert!(result.is_err());
}

#[test]
fn test_error_invalid_bond_order() {
    let sdf = r#"test


  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  9  0  0  0  0
M  END
$$$$
"#;
    // Bond order 9 is invalid
    let result = parse_sdf_string(sdf);
    assert!(result.is_err());
}
