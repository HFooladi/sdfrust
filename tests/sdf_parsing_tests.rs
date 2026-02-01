use sdfrust::{
    Atom, Bond, BondOrder, BondStereo, FileFormat, Molecule, detect_format, iter_auto_file,
    parse_auto_file, parse_auto_file_multi, parse_auto_string, parse_auto_string_multi,
    parse_sdf_file, parse_sdf_string, parse_sdf_string_multi, write_sdf_string,
};

const METHANE_SDF: &str = r#"methane
  sdfrust 3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;

const ETHANOL_SDF: &str = r#"ethanol
  sdfrust 3D

  9  8  0  0  0  0  0  0  0  0999 V2000
   -0.0014    1.0859    0.0082 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0018   -0.4224    0.0018 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2097   -0.9072   -0.0027 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.0197    1.4648    0.0014 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5363    1.4485   -0.8734 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5233    1.4374    0.9025 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5239   -0.7850    0.8911 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5103   -0.7927   -0.8837 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.1920   -1.8703   -0.0031 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
  1  6  1  0  0  0  0
  2  7  1  0  0  0  0
  2  8  1  0  0  0  0
  3  9  1  0  0  0  0
M  END
$$$$
"#;

const BENZENE_SDF: &str = r#"benzene
  sdfrust 3D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124   -0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2124    0.7000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  6  1  4  0  0  0  0
M  END
$$$$
"#;

const CHARGED_MOL_SDF: &str = r#"acetate
  sdfrust 3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.8000   -1.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  2  0  0  0  0
  2  4  1  0  0  0  0
M  CHG  1   4  -1
M  END
$$$$
"#;

// ============================================================================
// Basic parsing tests
// ============================================================================

#[test]
fn test_parse_methane() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_methane_atoms() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();

    // Carbon at origin
    let carbon = &mol.atoms[0];
    assert_eq!(carbon.element, "C");
    assert_eq!(carbon.x, 0.0);
    assert_eq!(carbon.y, 0.0);
    assert_eq!(carbon.z, 0.0);

    // Hydrogens
    for i in 1..5 {
        assert_eq!(mol.atoms[i].element, "H");
    }
}

#[test]
fn test_methane_bonds() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();

    // All bonds should be single bonds from carbon (index 0)
    for bond in &mol.bonds {
        assert_eq!(bond.order, BondOrder::Single);
        assert!(bond.atom1 == 0 || bond.atom2 == 0);
    }
}

#[test]
fn test_parse_ethanol() {
    let mol = parse_sdf_string(ETHANOL_SDF).unwrap();

    assert_eq!(mol.name, "ethanol");
    assert_eq!(mol.atom_count(), 9);
    assert_eq!(mol.bond_count(), 8);

    // Check element distribution
    let counts = mol.element_counts();
    assert_eq!(counts.get("C"), Some(&2));
    assert_eq!(counts.get("O"), Some(&1));
    assert_eq!(counts.get("H"), Some(&6));
}

#[test]
fn test_parse_benzene_aromatic() {
    let mol = parse_sdf_string(BENZENE_SDF).unwrap();

    assert_eq!(mol.name, "benzene");
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 6);

    // All bonds should be aromatic
    for bond in &mol.bonds {
        assert_eq!(bond.order, BondOrder::Aromatic);
    }

    assert!(mol.has_aromatic_bonds());
}

#[test]
fn test_parse_charged_molecule() {
    let mol = parse_sdf_string(CHARGED_MOL_SDF).unwrap();

    assert_eq!(mol.name, "acetate");

    // Check the charged oxygen
    let charged_o = &mol.atoms[3];
    assert_eq!(charged_o.element, "O");
    assert_eq!(charged_o.formal_charge, -1);

    // Total charge should be -1
    assert_eq!(mol.total_charge(), -1);
    assert!(mol.has_charges());
}

// ============================================================================
// Formula tests
// ============================================================================

#[test]
fn test_formula_methane() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();
    assert_eq!(mol.formula(), "CH4");
}

#[test]
fn test_formula_ethanol() {
    let mol = parse_sdf_string(ETHANOL_SDF).unwrap();
    assert_eq!(mol.formula(), "C2H6O");
}

#[test]
fn test_formula_benzene() {
    let mol = parse_sdf_string(BENZENE_SDF).unwrap();
    assert_eq!(mol.formula(), "C6");
}

// ============================================================================
// Geometry tests
// ============================================================================

#[test]
fn test_centroid() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();
    let (cx, cy, cz) = mol.centroid().unwrap();

    // Methane is symmetric around origin, centroid should be near (0, 0, 0)
    assert!(cx.abs() < 0.01);
    assert!(cy.abs() < 0.01);
    assert!(cz.abs() < 0.01);
}

#[test]
fn test_center_molecule() {
    let mut mol = parse_sdf_string(ETHANOL_SDF).unwrap();

    // Center the molecule
    mol.center();

    let (cx, cy, cz) = mol.centroid().unwrap();
    assert!(cx.abs() < 0.0001);
    assert!(cy.abs() < 0.0001);
    assert!(cz.abs() < 0.0001);
}

#[test]
fn test_translate() {
    let mut mol = parse_sdf_string(METHANE_SDF).unwrap();
    let original_x = mol.atoms[0].x;

    mol.translate(1.0, 2.0, 3.0);

    assert!((mol.atoms[0].x - (original_x + 1.0)).abs() < 0.0001);
    assert!((mol.atoms[0].y - 2.0).abs() < 0.0001);
    assert!((mol.atoms[0].z - 3.0).abs() < 0.0001);
}

#[test]
fn test_atom_distance() {
    let mol = parse_sdf_string(ETHANOL_SDF).unwrap();

    // C-C bond length should be around 1.5 Angstroms
    let c1 = &mol.atoms[0];
    let c2 = &mol.atoms[1];
    let cc_distance = c1.distance_to(c2);
    assert!(cc_distance > 1.4 && cc_distance < 1.6);

    // C-O bond length should be around 1.4 Angstroms
    let o = &mol.atoms[2];
    let co_distance = c2.distance_to(o);
    assert!(co_distance > 1.3 && co_distance < 1.5);
}

// ============================================================================
// Connectivity tests
// ============================================================================

#[test]
fn test_neighbors() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();

    // Carbon (index 0) should have 4 neighbors
    let carbon_neighbors = mol.neighbors(0);
    assert_eq!(carbon_neighbors.len(), 4);

    // Each hydrogen should have 1 neighbor (the carbon)
    for i in 1..5 {
        let h_neighbors = mol.neighbors(i);
        assert_eq!(h_neighbors.len(), 1);
        assert_eq!(h_neighbors[0], 0);
    }
}

#[test]
fn test_bonds_for_atom() {
    let mol = parse_sdf_string(BENZENE_SDF).unwrap();

    // Each carbon in benzene should have exactly 2 bonds (to neighboring carbons)
    for i in 0..6 {
        let bonds = mol.bonds_for_atom(i);
        assert_eq!(bonds.len(), 2);
    }
}

#[test]
fn test_bond_contains_atom() {
    let bond = Bond::new(0, 1, BondOrder::Single);
    assert!(bond.contains_atom(0));
    assert!(bond.contains_atom(1));
    assert!(!bond.contains_atom(2));
}

#[test]
fn test_bond_other_atom() {
    let bond = Bond::new(0, 1, BondOrder::Single);
    assert_eq!(bond.other_atom(0), Some(1));
    assert_eq!(bond.other_atom(1), Some(0));
    assert_eq!(bond.other_atom(2), None);
}

// ============================================================================
// Multi-molecule parsing tests
// ============================================================================

#[test]
fn test_parse_multi_molecule() {
    let multi_sdf = format!("{}{}{}", METHANE_SDF, ETHANOL_SDF, BENZENE_SDF);
    let mols = parse_sdf_string_multi(&multi_sdf).unwrap();

    assert_eq!(mols.len(), 3);
    assert_eq!(mols[0].name, "methane");
    assert_eq!(mols[1].name, "ethanol");
    assert_eq!(mols[2].name, "benzene");
}

// ============================================================================
// Properties tests
// ============================================================================

#[test]
fn test_parse_with_properties() {
    let sdf_with_props = r#"caffeine
  sdfrust 3D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
> <MOLECULAR_WEIGHT>
194.19

> <FORMULA>
C8H10N4O2

> <SMILES>
CN1C=NC2=C1C(=O)N(C(=O)N2C)C

$$$$
"#;

    let mol = parse_sdf_string(sdf_with_props).unwrap();

    assert_eq!(mol.get_property("MOLECULAR_WEIGHT"), Some("194.19"));
    assert_eq!(mol.get_property("FORMULA"), Some("C8H10N4O2"));
    assert_eq!(
        mol.get_property("SMILES"),
        Some("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    );
    assert_eq!(mol.get_property("NONEXISTENT"), None);
}

#[test]
fn test_set_property() {
    let mut mol = Molecule::new("test");
    mol.set_property("KEY1", "value1");
    mol.set_property("KEY2", "value2");

    assert_eq!(mol.get_property("KEY1"), Some("value1"));
    assert_eq!(mol.get_property("KEY2"), Some("value2"));

    // Overwrite
    mol.set_property("KEY1", "new_value");
    assert_eq!(mol.get_property("KEY1"), Some("new_value"));
}

// ============================================================================
// Writer tests
// ============================================================================

#[test]
fn test_write_and_parse_roundtrip() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();
    let output = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&output).unwrap();

    assert_eq!(parsed.name, mol.name);
    assert_eq!(parsed.atom_count(), mol.atom_count());
    assert_eq!(parsed.bond_count(), mol.bond_count());

    // Check coordinates
    for i in 0..mol.atom_count() {
        assert!((parsed.atoms[i].x - mol.atoms[i].x).abs() < 0.001);
        assert!((parsed.atoms[i].y - mol.atoms[i].y).abs() < 0.001);
        assert!((parsed.atoms[i].z - mol.atoms[i].z).abs() < 0.001);
    }
}

#[test]
fn test_write_preserves_properties() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.set_property("PROP1", "value1");
    mol.set_property("PROP2", "value2");

    let output = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&output).unwrap();

    assert_eq!(parsed.get_property("PROP1"), Some("value1"));
    assert_eq!(parsed.get_property("PROP2"), Some("value2"));
}

#[test]
fn test_write_aromatic_bonds() {
    let mol = parse_sdf_string(BENZENE_SDF).unwrap();
    let output = write_sdf_string(&mol).unwrap();
    let parsed = parse_sdf_string(&output).unwrap();

    // All bonds should still be aromatic
    for bond in &parsed.bonds {
        assert_eq!(bond.order, BondOrder::Aromatic);
    }
}

// ============================================================================
// Bond order tests
// ============================================================================

#[test]
fn test_bond_order_from_sdf() {
    assert_eq!(BondOrder::from_sdf(1), Some(BondOrder::Single));
    assert_eq!(BondOrder::from_sdf(2), Some(BondOrder::Double));
    assert_eq!(BondOrder::from_sdf(3), Some(BondOrder::Triple));
    assert_eq!(BondOrder::from_sdf(4), Some(BondOrder::Aromatic));
    assert_eq!(BondOrder::from_sdf(9), Some(BondOrder::Coordination));
    assert_eq!(BondOrder::from_sdf(10), Some(BondOrder::Hydrogen));
    assert_eq!(BondOrder::from_sdf(0), None);
    assert_eq!(BondOrder::from_sdf(11), None);
}

#[test]
fn test_bond_order_numeric() {
    assert_eq!(BondOrder::Single.order(), 1.0);
    assert_eq!(BondOrder::Double.order(), 2.0);
    assert_eq!(BondOrder::Triple.order(), 3.0);
    assert_eq!(BondOrder::Aromatic.order(), 1.5);
}

#[test]
fn test_bond_stereo() {
    assert_eq!(BondStereo::from_sdf(0), BondStereo::None);
    assert_eq!(BondStereo::from_sdf(1), BondStereo::Up);
    assert_eq!(BondStereo::from_sdf(4), BondStereo::Either);
    assert_eq!(BondStereo::from_sdf(6), BondStereo::Down);
}

// ============================================================================
// Edge cases and error handling
// ============================================================================

#[test]
fn test_empty_molecule_name() {
    let sdf = r#"
  sdfrust 3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "");
    assert_eq!(mol.atom_count(), 1);
}

#[test]
fn test_molecule_with_no_bonds() {
    let sdf = r#"helium
  sdfrust 3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 He  0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

    let mol = parse_sdf_string(sdf).unwrap();
    assert_eq!(mol.name, "helium");
    assert_eq!(mol.atom_count(), 1);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(mol.formula(), "He");
}

#[test]
fn test_total_bond_order() {
    let mol = parse_sdf_string(METHANE_SDF).unwrap();
    assert_eq!(mol.total_bond_order(), 4.0); // 4 single bonds

    let benzene = parse_sdf_string(BENZENE_SDF).unwrap();
    assert_eq!(benzene.total_bond_order(), 9.0); // 6 aromatic bonds (1.5 each)
}

#[test]
fn test_atoms_by_element() {
    let mol = parse_sdf_string(ETHANOL_SDF).unwrap();

    let carbons = mol.atoms_by_element("C");
    assert_eq!(carbons.len(), 2);

    let oxygens = mol.atoms_by_element("O");
    assert_eq!(oxygens.len(), 1);

    let hydrogens = mol.atoms_by_element("H");
    assert_eq!(hydrogens.len(), 6);
}

#[test]
fn test_bonds_by_order() {
    let mol = parse_sdf_string(CHARGED_MOL_SDF).unwrap();

    let single_bonds = mol.bonds_by_order(BondOrder::Single);
    assert_eq!(single_bonds.len(), 2);

    let double_bonds = mol.bonds_by_order(BondOrder::Double);
    assert_eq!(double_bonds.len(), 1);
}

// ============================================================================
// Real-world file tests
// ============================================================================

#[test]
fn test_parse_caffeine_from_file() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();

    // PubChem uses CID as molecule name
    assert_eq!(mol.name, "2519");
    assert_eq!(mol.atom_count(), 24);
    assert_eq!(mol.bond_count(), 25);

    // Check element counts (C8H10N4O2)
    let counts = mol.element_counts();
    assert_eq!(counts.get("C"), Some(&8));
    assert_eq!(counts.get("H"), Some(&10));
    assert_eq!(counts.get("N"), Some(&4));
    assert_eq!(counts.get("O"), Some(&2));

    // Check PubChem properties were parsed
    assert_eq!(mol.get_property("PUBCHEM_COMPOUND_CID"), Some("2519"));
    assert_eq!(
        mol.get_property("PUBCHEM_MOLECULAR_FORMULA"),
        Some("C8H10N4O2")
    );
    assert_eq!(mol.get_property("PUBCHEM_MOLECULAR_WEIGHT"), Some("194.19"));

    // Verify our formula calculation matches PubChem
    assert_eq!(mol.formula(), "C8H10N4O2");
}

// ============================================================================
// Auto-detection tests (FileFormat, detect_format, parse_auto_*)
// ============================================================================

const MOL2_CONTENT: &str = r#"@<TRIPOS>MOLECULE
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

const V3000_CONTENT: &str = r#"methane
  sdfrust 3D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 H 0.6289 0.6289 0.6289 0
M  V30 3 H -0.6289 -0.6289 0.6289 0
M  V30 4 H -0.6289 0.6289 -0.6289 0
M  V30 5 H 0.6289 -0.6289 -0.6289 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;

#[test]
fn test_detect_format_mol2() {
    assert_eq!(detect_format(MOL2_CONTENT), FileFormat::Mol2);
}

#[test]
fn test_detect_format_v3000() {
    assert_eq!(detect_format(V3000_CONTENT), FileFormat::SdfV3000);
}

#[test]
fn test_detect_format_v2000() {
    assert_eq!(detect_format(METHANE_SDF), FileFormat::SdfV2000);
    assert_eq!(detect_format(ETHANOL_SDF), FileFormat::SdfV2000);
    assert_eq!(detect_format(BENZENE_SDF), FileFormat::SdfV2000);
}

#[test]
fn test_detect_format_mol2_marker_deep() {
    // Test that we detect MOL2 even if @<TRIPOS> appears after some lines
    let content = "some header\nanother line\n@<TRIPOS>MOLECULE\ntest\n";
    assert_eq!(detect_format(content), FileFormat::Mol2);
}

#[test]
fn test_detect_format_empty_defaults_v2000() {
    // Empty or minimal content should default to V2000
    assert_eq!(detect_format(""), FileFormat::SdfV2000);
    assert_eq!(detect_format("test"), FileFormat::SdfV2000);
}

#[test]
fn test_file_format_display() {
    assert_eq!(format!("{}", FileFormat::SdfV2000), "sdf_v2000");
    assert_eq!(format!("{}", FileFormat::SdfV3000), "sdf_v3000");
    assert_eq!(format!("{}", FileFormat::Mol2), "mol2");
}

#[test]
fn test_parse_auto_string_v2000() {
    let mol = parse_auto_string(METHANE_SDF).unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_parse_auto_string_v3000() {
    let mol = parse_auto_string(V3000_CONTENT).unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_parse_auto_string_mol2() {
    let mol = parse_auto_string(MOL2_CONTENT).unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_parse_auto_string_multi_v2000() {
    let multi_sdf = format!("{}{}", METHANE_SDF, ETHANOL_SDF);
    let mols = parse_auto_string_multi(&multi_sdf).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "methane");
    assert_eq!(mols[1].name, "ethanol");
}

#[test]
fn test_parse_auto_string_multi_mol2() {
    let multi_mol2 = format!("{}{}", MOL2_CONTENT, MOL2_CONTENT);
    let mols = parse_auto_string_multi(&multi_mol2).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "methane");
    assert_eq!(mols[1].name, "methane");
}

#[test]
fn test_parse_auto_file_v2000() {
    let mol = parse_auto_file("tests/test_data/aspirin.sdf").unwrap();
    assert_eq!(mol.name, "2244");
    assert_eq!(mol.atom_count(), 21);
}

#[test]
fn test_parse_auto_file_v3000() {
    let mol = parse_auto_file("tests/test_data/v3000_methane.sdf").unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
}

#[test]
fn test_parse_auto_file_mol2() {
    let mol = parse_auto_file("tests/test_data/methane.mol2").unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
}

#[test]
fn test_parse_auto_file_multi_v2000() {
    let mols = parse_auto_file_multi("tests/test_data/caffeine_pubchem.sdf").unwrap();
    assert!(!mols.is_empty());
    assert_eq!(mols[0].name, "2519");
}

#[test]
fn test_iter_auto_file_v2000() {
    let iter = iter_auto_file("tests/test_data/aspirin.sdf").unwrap();
    let mols: Vec<_> = iter.collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].name, "2244");
}

#[test]
fn test_iter_auto_file_mol2() {
    let iter = iter_auto_file("tests/test_data/methane.mol2").unwrap();
    let mols: Vec<_> = iter.collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(mols.len(), 1);
    assert_eq!(mols[0].name, "methane");
}

#[test]
fn test_iter_auto_file_v3000() {
    let iter = iter_auto_file("tests/test_data/v3000_benzene.sdf").unwrap();
    let mols: Vec<_> = iter.collect::<Result<Vec<_>, _>>().unwrap();
    assert!(!mols.is_empty());
    assert_eq!(mols[0].name, "benzene");
}
