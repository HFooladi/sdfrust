//! Tests for SDF V3000 format parsing and writing.

use sdfrust::{
    Atom, Bond, BondOrder, BondStereo, Molecule, SdfFormat, StereoGroup, detect_sdf_format,
    needs_v3000, parse_sdf_auto_string, parse_sdf_v3000_file, parse_sdf_v3000_string,
    write_sdf_auto_string, write_sdf_v3000_string,
};

// ============================================================================
// Basic V3000 Parsing Tests
// ============================================================================

#[test]
fn test_parse_v3000_methane_file() {
    let mol = parse_sdf_v3000_file("tests/test_data/v3000_methane.sdf").unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
    assert_eq!(mol.formula(), "CH4");
    assert_eq!(mol.format_version, SdfFormat::V3000);

    // Check properties
    assert_eq!(mol.get_property("FORMULA"), Some("CH4"));
    assert_eq!(mol.get_property("MW"), Some("16.04"));
}

#[test]
fn test_parse_v3000_benzene_file() {
    let mol = parse_sdf_v3000_file("tests/test_data/v3000_benzene.sdf").unwrap();

    assert_eq!(mol.name, "benzene");
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);
    assert_eq!(mol.formula(), "C6");
    assert!(mol.has_aromatic_bonds());

    // All bonds should be aromatic
    for bond in &mol.bonds {
        assert_eq!(bond.order, BondOrder::Aromatic);
    }
}

#[test]
fn test_parse_v3000_charged_file() {
    let mol = parse_sdf_v3000_file("tests/test_data/v3000_charged.sdf").unwrap();

    assert_eq!(mol.name, "ammonium_acetate");
    assert!(mol.has_charges());

    // Find the nitrogen (should have +1 charge)
    let nitrogen = mol.atoms.iter().find(|a| a.element == "N").unwrap();
    assert_eq!(nitrogen.formal_charge, 1);

    // Find one of the oxygens (should have -1 charge)
    let charged_oxygen = mol
        .atoms
        .iter()
        .find(|a| a.element == "O" && a.formal_charge == -1);
    assert!(charged_oxygen.is_some());

    // Total charge should be 0
    assert_eq!(mol.total_charge(), 0);
}

#[test]
fn test_parse_v3000_stereo_file() {
    let mol = parse_sdf_v3000_file("tests/test_data/v3000_stereo.sdf").unwrap();

    assert_eq!(mol.name, "chiral_molecule");
    assert_eq!(mol.atom_count(), 5);

    // Check stereo bonds
    let wedge_bond = mol.bonds.iter().find(|b| b.stereo == BondStereo::Up);
    assert!(wedge_bond.is_some());

    let dash_bond = mol.bonds.iter().find(|b| b.stereo == BondStereo::Down);
    assert!(dash_bond.is_some());

    // Check property
    assert_eq!(mol.get_property("STEREOCHEMISTRY"), Some("absolute"));
}

// ============================================================================
// V3000 String Parsing Tests
// ============================================================================

#[test]
fn test_parse_v3000_string() {
    let content = r#"water
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0000 0.0000 0.0000 0
M  V30 2 H 0.9572 0.0000 0.0000 0
M  V30 3 H -0.2400 0.9266 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.name, "water");
    assert_eq!(mol.formula(), "H2O");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_parse_v3000_with_v3000_ids() {
    let content = r#"test
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 10 C 0.0000 0.0000 0.0000 0
M  V30 20 O 1.2000 0.0000 0.0000 0
M  V30 30 O -1.2000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 10 20
M  V30 2 2 10 30
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    // Check that non-sequential IDs are handled correctly
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);

    // Atoms should have their V3000 IDs stored
    assert_eq!(mol.atoms[0].v3000_id, Some(10));
    assert_eq!(mol.atoms[1].v3000_id, Some(20));
    assert_eq!(mol.atoms[2].v3000_id, Some(30));

    // Bonds should reference correct atoms by index
    assert_eq!(mol.bonds[0].atom1, 0); // atom with V3000 ID 10
    assert_eq!(mol.bonds[0].atom2, 1); // atom with V3000 ID 20
}

#[test]
fn test_parse_v3000_double_bond() {
    let content = r#"formaldehyde
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 O 1.2000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.bonds[0].order, BondOrder::Double);
}

#[test]
fn test_parse_v3000_triple_bond() {
    let content = r#"hydrogen_cyanide
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 H 0.0000 0.0000 0.0000 0
M  V30 2 C 1.0000 0.0000 0.0000 0
M  V30 3 N 2.1500 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 3 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    let triple_bond = mol.bonds.iter().find(|b| b.order == BondOrder::Triple);
    assert!(triple_bond.is_some());
}

// ============================================================================
// V3000 Extended Bond Types Tests
// ============================================================================

#[test]
fn test_parse_v3000_coordination_bond() {
    let content = r#"coordination
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0000 0.0000 0.0000 0
M  V30 2 Fe 2.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 9 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.bonds[0].order, BondOrder::Coordination);
}

#[test]
fn test_parse_v3000_hydrogen_bond() {
    let content = r#"hbond
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O 0.0000 0.0000 0.0000 0
M  V30 2 H 0.9572 0.0000 0.0000 0
M  V30 3 O 2.5000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 10 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    let hbond = mol.bonds.iter().find(|b| b.order == BondOrder::Hydrogen);
    assert!(hbond.is_some());
}

// ============================================================================
// V3000 Charge and Radical Tests
// ============================================================================

#[test]
fn test_parse_v3000_charges() {
    let content = r#"charged
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0000 0.0000 0.0000 0 CHG=1
M  V30 2 O 1.5000 0.0000 0.0000 0 CHG=-1
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.atoms[0].formal_charge, 1);
    assert_eq!(mol.atoms[1].formal_charge, -1);
    assert_eq!(mol.total_charge(), 0);
}

#[test]
fn test_parse_v3000_radical() {
    let content = r#"radical
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 2 1 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0 RAD=2
M  V30 2 H 1.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.atoms[0].radical, Some(2)); // doublet
}

// ============================================================================
// Format Detection Tests
// ============================================================================

#[test]
fn test_detect_v2000_format() {
    let content = r#"test
  sdfrust   01012500003D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;
    assert_eq!(detect_sdf_format(content), SdfFormat::V2000);
}

#[test]
fn test_detect_v3000_format() {
    let content = r#"test
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
"#;
    assert_eq!(detect_sdf_format(content), SdfFormat::V3000);
}

#[test]
fn test_parse_auto_v2000() {
    let content = r#"methane
  sdfrust   01012500003D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;
    let mol = parse_sdf_auto_string(content).unwrap();
    assert_eq!(mol.format_version, SdfFormat::V2000);
}

#[test]
fn test_parse_auto_v3000() {
    let content = r#"methane
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_auto_string(content).unwrap();
    assert_eq!(mol.format_version, SdfFormat::V3000);
}

// ============================================================================
// V3000 Writing Tests
// ============================================================================

#[test]
fn test_write_v3000_simple() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

    let output = write_sdf_v3000_string(&mol).unwrap();

    assert!(output.contains("V3000"));
    assert!(output.contains("BEGIN CTAB"));
    assert!(output.contains("COUNTS 2 1"));
    assert!(output.contains("BEGIN ATOM"));
    assert!(output.contains("END ATOM"));
    assert!(output.contains("BEGIN BOND"));
    assert!(output.contains("END BOND"));
    assert!(output.contains("END CTAB"));
    assert!(output.contains("M  END"));
    assert!(output.contains("$$$$"));
}

#[test]
fn test_write_v3000_with_charge() {
    let mut mol = Molecule::new("charged");
    mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
    mol.atoms[0].formal_charge = 1;
    mol.atoms.push(Atom::new(1, "Cl", 2.0, 0.0, 0.0));
    mol.atoms[1].formal_charge = -1;
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

    let output = write_sdf_v3000_string(&mol).unwrap();

    assert!(output.contains("CHG=1"));
    assert!(output.contains("CHG=-1"));
}

#[test]
fn test_write_v3000_with_properties() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.set_property("SMILES", "C");
    mol.set_property("MW", "12.01");

    let output = write_sdf_v3000_string(&mol).unwrap();

    assert!(output.contains("> <SMILES>"));
    assert!(output.contains("C"));
    assert!(output.contains("> <MW>"));
    assert!(output.contains("12.01"));
}

// ============================================================================
// Round-Trip Tests
// ============================================================================

#[test]
fn test_v3000_round_trip() {
    let mut mol = Molecule::new("round_trip");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "O", 1.2, 0.5, 0.0));
    mol.atoms.push(Atom::new(2, "O", -1.2, -0.5, 0.0));
    mol.atoms[1].formal_charge = -1;
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 2, BondOrder::Double));
    mol.set_property("NAME", "formate");

    let v3000_str = write_sdf_v3000_string(&mol).unwrap();
    let parsed = parse_sdf_v3000_string(&v3000_str).unwrap();

    assert_eq!(parsed.name, mol.name);
    assert_eq!(parsed.atom_count(), mol.atom_count());
    assert_eq!(parsed.bond_count(), mol.bond_count());
    assert_eq!(parsed.atoms[1].formal_charge, -1);
    assert_eq!(parsed.get_property("NAME"), Some("formate"));
}

#[test]
fn test_v3000_round_trip_aromatics() {
    let content = std::fs::read_to_string("tests/test_data/v3000_benzene.sdf").unwrap();
    let mol = parse_sdf_v3000_string(&content).unwrap();

    let v3000_str = write_sdf_v3000_string(&mol).unwrap();
    let parsed = parse_sdf_v3000_string(&v3000_str).unwrap();

    assert_eq!(parsed.name, mol.name);
    assert_eq!(parsed.atom_count(), mol.atom_count());
    assert_eq!(parsed.bond_count(), mol.bond_count());
    assert!(parsed.has_aromatic_bonds());
}

// ============================================================================
// Needs V3000 Detection Tests
// ============================================================================

#[test]
fn test_needs_v3000_small_molecule() {
    let mut mol = Molecule::new("small");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 0, BondOrder::Single));

    assert!(!needs_v3000(&mol));
}

#[test]
fn test_needs_v3000_coordination_bond() {
    let mut mol = Molecule::new("coordination");
    mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "Fe", 2.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Coordination));

    assert!(needs_v3000(&mol));
}

#[test]
fn test_needs_v3000_stereogroups() {
    let mut mol = Molecule::new("stereo");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.stereogroups.push(StereoGroup::absolute(vec![0]));

    assert!(needs_v3000(&mol));
}

// ============================================================================
// Auto Format Selection Tests
// ============================================================================

#[test]
fn test_write_auto_uses_v2000_for_simple() {
    let mut mol = Molecule::new("simple");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 0, BondOrder::Single));

    let output = write_sdf_auto_string(&mol).unwrap();

    assert!(output.contains("V2000"));
    assert!(!output.contains("V3000"));
}

#[test]
fn test_write_auto_uses_v3000_for_extended() {
    let mut mol = Molecule::new("extended");
    mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "Fe", 2.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Coordination));

    let output = write_sdf_auto_string(&mol).unwrap();

    assert!(output.contains("V3000"));
    assert!(!output.contains("V2000"));
}

#[test]
fn test_write_auto_uses_v3000_for_v3000_format() {
    let mut mol = Molecule::new("v3000_format");
    mol.format_version = SdfFormat::V3000;
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

    let output = write_sdf_auto_string(&mol).unwrap();

    assert!(output.contains("V3000"));
}

// ============================================================================
// Edge Cases and Error Handling
// ============================================================================

#[test]
fn test_parse_v3000_empty_molecule() {
    let content = r#"empty
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 0 0 0 0 0
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.name, "empty");
    assert_eq!(mol.atom_count(), 0);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_parse_v3000_atoms_only() {
    let content = r#"atoms_only
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 He 0.0000 0.0000 0.0000 0
M  V30 2 Ne 2.0000 0.0000 0.0000 0
M  V30 3 Ar 4.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
"#;
    let mol = parse_sdf_v3000_string(content).unwrap();

    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_v3000_coordinate_precision() {
    let mut mol = Molecule::new("precision");
    mol.atoms
        .push(Atom::new(0, "C", 1.23456789, -2.34567891, 3.45678912));

    let output = write_sdf_v3000_string(&mol).unwrap();
    let parsed = parse_sdf_v3000_string(&output).unwrap();

    // Check coordinates are preserved to 4 decimal places
    let atom = &parsed.atoms[0];
    assert!((atom.x - 1.2346).abs() < 0.0001);
    assert!((atom.y - (-2.3457)).abs() < 0.0001);
    assert!((atom.z - 3.4568).abs() < 0.0001);
}

// ============================================================================
// Multi-molecule Tests
// ============================================================================

#[test]
fn test_parse_v3000_multi() {
    let content = r#"mol1
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
mol2
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N 0.0000 0.0000 0.0000 0
M  V30 END ATOM
M  V30 END CTAB
M  END
$$$$
"#;
    let mols = sdfrust::parse_sdf_v3000_string_multi(content).unwrap();

    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "mol1");
    assert_eq!(mols[0].atoms[0].element, "C");
    assert_eq!(mols[1].name, "mol2");
    assert_eq!(mols[1].atoms[0].element, "N");
}
