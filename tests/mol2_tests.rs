//! Integration tests for MOL2 parsing.
//!
//! These tests verify the MOL2 parser against real file formats.

use sdfrust::{parse_mol2_file, parse_mol2_string, parse_mol2_string_multi, BondOrder};

const TEST_DATA_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test_data");

// ============================================================================
// File Parsing Tests
// ============================================================================

#[test]
fn test_parse_methane_mol2_file() {
    let path = format!("{}/methane.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
    assert_eq!(mol.formula(), "CH4");
}

#[test]
fn test_parse_benzene_mol2_file() {
    let path = format!("{}/benzene.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    assert_eq!(mol.name, "benzene");
    assert_eq!(mol.atom_count(), 12); // 6 C + 6 H
    assert_eq!(mol.bond_count(), 12); // 6 aromatic + 6 C-H
    assert_eq!(mol.formula(), "C6H6");
}

#[test]
fn test_benzene_aromatic_bonds() {
    let path = format!("{}/benzene.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    // Count aromatic vs single bonds
    let aromatic_count = mol.bonds.iter().filter(|b| b.order == BondOrder::Aromatic).count();
    let single_count = mol.bonds.iter().filter(|b| b.order == BondOrder::Single).count();

    assert_eq!(aromatic_count, 6, "Benzene should have 6 aromatic bonds");
    assert_eq!(single_count, 6, "Benzene should have 6 single C-H bonds");
}

#[test]
fn test_benzene_carbon_positions() {
    let path = format!("{}/benzene.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    // All carbons should be in z=0 plane
    let carbons: Vec<_> = mol.atoms.iter().filter(|a| a.element == "C").collect();
    assert_eq!(carbons.len(), 6);

    for c in carbons {
        assert!((c.z - 0.0).abs() < 1e-6, "Carbon z coordinate should be 0");
    }
}

// ============================================================================
// String Parsing Tests
// ============================================================================

#[test]
fn test_parse_water_mol2_string() {
    let water_mol2 = r#"@<TRIPOS>MOLECULE
water
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
"#;

    let mol = parse_mol2_string(water_mol2).unwrap();

    assert_eq!(mol.name, "water");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
    assert_eq!(mol.formula(), "H2O");

    // Check oxygen is at origin
    assert_eq!(mol.atoms[0].element, "O");
    assert_eq!(mol.atoms[0].x, 0.0);
    assert_eq!(mol.atoms[0].y, 0.0);
    assert_eq!(mol.atoms[0].z, 0.0);
}

#[test]
fn test_parse_ethane_mol2_string() {
    let ethane_mol2 = r#"@<TRIPOS>MOLECULE
ethane
 8 7 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 C2          1.5400    0.0000    0.0000 C.3       1 MOL       0.0000
      3 H1         -0.3700    1.0260    0.0000 H         1 MOL       0.0000
      4 H2         -0.3700   -0.5130    0.8890 H         1 MOL       0.0000
      5 H3         -0.3700   -0.5130   -0.8890 H         1 MOL       0.0000
      6 H4          1.9100    1.0260    0.0000 H         1 MOL       0.0000
      7 H5          1.9100   -0.5130    0.8890 H         1 MOL       0.0000
      8 H6          1.9100   -0.5130   -0.8890 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
     5     2     6 1
     6     2     7 1
     7     2     8 1
"#;

    let mol = parse_mol2_string(ethane_mol2).unwrap();

    assert_eq!(mol.name, "ethane");
    assert_eq!(mol.atom_count(), 8);
    assert_eq!(mol.bond_count(), 7);
    assert_eq!(mol.formula(), "C2H6");

    // All bonds should be single
    for bond in &mol.bonds {
        assert_eq!(bond.order, BondOrder::Single);
    }
}

#[test]
fn test_parse_double_bond_mol2() {
    let ethene_mol2 = r#"@<TRIPOS>MOLECULE
ethene
 6 5 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.2       1 MOL       0.0000
      2 C2          1.3300    0.0000    0.0000 C.2       1 MOL       0.0000
      3 H1         -0.5200    0.9430    0.0000 H         1 MOL       0.0000
      4 H2         -0.5200   -0.9430    0.0000 H         1 MOL       0.0000
      5 H3          1.8500    0.9430    0.0000 H         1 MOL       0.0000
      6 H4          1.8500   -0.9430    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 2
     2     1     3 1
     3     1     4 1
     4     2     5 1
     5     2     6 1
"#;

    let mol = parse_mol2_string(ethene_mol2).unwrap();

    assert_eq!(mol.name, "ethene");
    assert_eq!(mol.formula(), "C2H4");

    // First bond should be double
    assert_eq!(mol.bonds[0].order, BondOrder::Double);

    // Rest should be single
    for bond in mol.bonds.iter().skip(1) {
        assert_eq!(bond.order, BondOrder::Single);
    }
}

#[test]
fn test_parse_triple_bond_mol2() {
    let ethyne_mol2 = r#"@<TRIPOS>MOLECULE
ethyne
 4 3 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.1       1 MOL       0.0000
      2 C2          1.2000    0.0000    0.0000 C.1       1 MOL       0.0000
      3 H1         -1.0600    0.0000    0.0000 H         1 MOL       0.0000
      4 H2          2.2600    0.0000    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 3
     2     1     3 1
     3     2     4 1
"#;

    let mol = parse_mol2_string(ethyne_mol2).unwrap();

    assert_eq!(mol.name, "ethyne");
    assert_eq!(mol.formula(), "C2H2");

    // First bond should be triple
    assert_eq!(mol.bonds[0].order, BondOrder::Triple);
}

// ============================================================================
// Multi-molecule Tests
// ============================================================================

#[test]
fn test_parse_multi_mol2() {
    let multi_mol2 = r#"@<TRIPOS>MOLECULE
water
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
@<TRIPOS>MOLECULE
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

    let mols = parse_mol2_string_multi(multi_mol2).unwrap();

    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "water");
    assert_eq!(mols[0].formula(), "H2O");
    assert_eq!(mols[1].name, "methane");
    assert_eq!(mols[1].formula(), "CH4");
}

// ============================================================================
// Charged Molecule Tests
// ============================================================================

#[test]
fn test_parse_charged_molecule_mol2() {
    let acetate_mol2 = r#"@<TRIPOS>MOLECULE
acetate
 4 3 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 C2          1.5000    0.0000    0.0000 C.2       1 MOL       0.8000
      3 O1          2.0500    1.0700    0.0000 O.co2     1 MOL      -0.9000
      4 O2          2.0500   -1.0700    0.0000 O.co2     1 MOL      -0.9000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 1
     3     2     4 1
"#;

    let mol = parse_mol2_string(acetate_mol2).unwrap();

    assert_eq!(mol.name, "acetate");
    assert_eq!(mol.atom_count(), 4);

    // Check that partial charges are converted to formal charges
    assert_eq!(mol.atoms[1].formal_charge, 1); // C2 has +0.8, rounds to +1
    assert_eq!(mol.atoms[2].formal_charge, -1); // O1 has -0.9, rounds to -1
    assert_eq!(mol.atoms[3].formal_charge, -1); // O2 has -0.9, rounds to -1
}

// ============================================================================
// Atom Type Tests
// ============================================================================

#[test]
fn test_mol2_atom_type_parsing() {
    let mol2 = r#"@<TRIPOS>MOLECULE
test
 3 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.ar      1 MOL       0.0000
      2 N1          1.0000    0.0000    0.0000 N.pl3     1 MOL       0.0000
      3 O1          2.0000    0.0000    0.0000 O.co2     1 MOL       0.0000
@<TRIPOS>BOND
"#;

    let mol = parse_mol2_string(mol2).unwrap();

    // Atom types like "C.ar" should be parsed as element "C"
    assert_eq!(mol.atoms[0].element, "C");
    assert_eq!(mol.atoms[1].element, "N");
    assert_eq!(mol.atoms[2].element, "O");
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn test_mol2_with_extra_sections() {
    let mol2_with_substructure = r#"@<TRIPOS>MOLECULE
test
 3 2 1 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 RESIDUE           4 A     MOL     0 ROOT
"#;

    let mol = parse_mol2_string(mol2_with_substructure).unwrap();

    // Should still parse correctly, ignoring SUBSTRUCTURE section
    assert_eq!(mol.name, "test");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_mol2_with_comment_section() {
    let mol2_with_comment = r#"@<TRIPOS>MOLECULE
test_molecule
 5 4 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>COMMENT
This is a test comment that spans
multiple lines and should be ignored.

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

    let mol = parse_mol2_string(mol2_with_comment).unwrap();

    assert_eq!(mol.name, "test_molecule");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

// ============================================================================
// Molecule Method Tests
// ============================================================================

#[test]
fn test_mol2_centroid_calculation() {
    let path = format!("{}/methane.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    let centroid = mol.centroid();

    // Carbon is at origin, hydrogens symmetrically distributed
    // Centroid should be at the center of mass
    assert!(centroid.is_some());
    let (cx, cy, cz) = centroid.unwrap();

    // For methane with C at origin, centroid should be very close to origin
    assert!((cx).abs() < 0.1, "Centroid x should be near 0");
    assert!((cy).abs() < 0.1, "Centroid y should be near 0");
    assert!((cz).abs() < 0.1, "Centroid z should be near 0");
}

#[test]
fn test_mol2_element_counts() {
    let path = format!("{}/benzene.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    let counts = mol.element_counts();

    assert_eq!(counts.get("C"), Some(&6));
    assert_eq!(counts.get("H"), Some(&6));
}

#[test]
fn test_mol2_neighbors() {
    let path = format!("{}/methane.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    // Carbon (index 0) should have 4 hydrogen neighbors
    let carbon_neighbors = mol.neighbors(0);
    assert_eq!(carbon_neighbors.len(), 4);

    // Each hydrogen should have 1 neighbor (the carbon)
    for i in 1..5 {
        let h_neighbors = mol.neighbors(i);
        assert_eq!(h_neighbors.len(), 1);
        assert_eq!(h_neighbors[0], 0); // Connected to carbon
    }
}

#[test]
fn test_mol2_total_bond_order() {
    let path = format!("{}/benzene.mol2", TEST_DATA_DIR);
    let mol = parse_mol2_file(&path).unwrap();

    // 6 aromatic bonds (1.5 each) + 6 single bonds = 6*1.5 + 6*1 = 15
    let total = mol.total_bond_order();
    assert!((total - 15.0).abs() < 0.01);
}
