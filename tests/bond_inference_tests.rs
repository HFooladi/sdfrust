//! Integration tests for bond inference from 3D coordinates.

use sdfrust::{
    Atom, Bond, BondInferenceConfig, BondOrder, Molecule, infer_bonds, infer_bonds_with_config,
    parse_xyz_file, parse_xyz_file_multi, parse_xyz_string,
};

#[test]
fn test_water_xyz_file() {
    let mut mol = parse_xyz_file("tests/test_data/water.xyz").unwrap();
    assert_eq!(mol.bond_count(), 0);

    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 2);

    // Both bonds should be O-H
    for bond in &mol.bonds {
        let elem1 = &mol.atoms[bond.atom1].element;
        let elem2 = &mol.atoms[bond.atom2].element;
        assert!(
            (elem1 == "O" && elem2 == "H") || (elem1 == "H" && elem2 == "O"),
            "Expected O-H bond, got {}-{}",
            elem1,
            elem2
        );
        assert_eq!(bond.order, BondOrder::Single);
    }
}

#[test]
fn test_multi_xyz_file() {
    let mut mols = parse_xyz_file_multi("tests/test_data/multi.xyz").unwrap();
    assert_eq!(mols.len(), 3);

    // Water: 2 bonds
    infer_bonds(&mut mols[0], None).unwrap();
    assert_eq!(mols[0].bond_count(), 2);

    // Methane: 4 bonds
    infer_bonds(&mut mols[1], None).unwrap();
    assert_eq!(mols[1].bond_count(), 4);

    // H2: 1 bond
    infer_bonds(&mut mols[2], None).unwrap();
    assert_eq!(mols[2].bond_count(), 1);
}

#[test]
fn test_ethanol_inline() {
    // Ethanol (C2H5OH) - 8 bonds: 1 C-C, 1 C-O, 5 C-H, 1 O-H
    let xyz = "9
ethanol
C   0.0000   0.0000   0.0000
C   1.5200   0.0000   0.0000
O   2.0800   1.2100   0.0000
H  -0.3600   1.0100   0.0000
H  -0.3600  -0.5100   0.8800
H  -0.3600  -0.5100  -0.8800
H   1.8800  -0.5100   0.8800
H   1.8800  -0.5100  -0.8800
H   3.0400   1.2100   0.0000
";
    let mut mol = parse_xyz_string(xyz).unwrap();
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 8);
}

#[test]
fn test_co2_inline() {
    // CO2 - 2 bonds (C=O, but inferred as single)
    let xyz = "3
carbon dioxide
C   0.0000   0.0000   0.0000
O   1.1600   0.0000   0.0000
O  -1.1600   0.0000   0.0000
";
    let mut mol = parse_xyz_string(xyz).unwrap();
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_tolerance_effect() {
    // Two C atoms at distance 1.9 A
    // Sum of covalent radii = 0.77 + 0.77 = 1.54
    // Default tolerance (0.45): 1.54 + 0.45 = 1.99 → bonded
    // Tight tolerance (0.1): 1.54 + 0.1 = 1.64 → not bonded
    let xyz = "2
test
C   0.0000   0.0000   0.0000
C   1.9000   0.0000   0.0000
";
    let mut mol = parse_xyz_string(xyz).unwrap();

    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(
        mol.bond_count(),
        1,
        "Should be bonded with default tolerance"
    );

    infer_bonds(&mut mol, Some(0.1)).unwrap();
    assert_eq!(
        mol.bond_count(),
        0,
        "Should not be bonded with tight tolerance"
    );
}

#[test]
fn test_empty_molecule() {
    let mut mol = Molecule::new("empty");
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_single_atom() {
    let mut mol = Molecule::new("single");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_unknown_element_error() {
    let mut mol = Molecule::new("unknown");
    mol.atoms.push(Atom::new(0, "Xx", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    let result = infer_bonds(&mut mol, None);
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("Xx"));
    assert!(err_msg.contains("index 0"));
}

#[test]
fn test_distant_atoms_no_bond() {
    let mut mol = Molecule::new("distant");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 100.0, 0.0, 0.0));
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 0);
}

#[test]
fn test_clears_existing_bonds() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Triple));

    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bonds[0].order, BondOrder::Single);
}

#[test]
fn test_config_keep_existing_bonds() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

    let config = BondInferenceConfig {
        clear_existing_bonds: false,
        ..Default::default()
    };
    infer_bonds_with_config(&mut mol, &config).unwrap();
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_molecule_method() {
    let mut mol = Molecule::new("water");
    mol.atoms
        .push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
    mol.atoms
        .push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
    mol.atoms
        .push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));

    mol.infer_bonds(None).unwrap();
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_molecule_method_with_tolerance() {
    let mut mol = Molecule::new("water");
    mol.atoms
        .push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
    mol.atoms
        .push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
    mol.atoms
        .push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));

    mol.infer_bonds(Some(0.5)).unwrap();
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_benzene_ring() {
    // Benzene with approximate 3D coordinates
    let xyz = "12
benzene
C   1.2124   0.7000   0.0000
C   1.2124  -0.7000   0.0000
C   0.0000  -1.4000   0.0000
C  -1.2124  -0.7000   0.0000
C  -1.2124   0.7000   0.0000
C   0.0000   1.4000   0.0000
H   2.1560   1.2450   0.0000
H   2.1560  -1.2450   0.0000
H   0.0000  -2.4900   0.0000
H  -2.1560  -1.2450   0.0000
H  -2.1560   1.2450   0.0000
H   0.0000   2.4900   0.0000
";
    let mut mol = parse_xyz_string(xyz).unwrap();
    infer_bonds(&mut mol, None).unwrap();
    // 6 C-C + 6 C-H = 12 bonds
    assert_eq!(mol.bond_count(), 12);
}

#[test]
fn test_overlapping_atoms_no_bond() {
    // Two atoms at the same position should not be bonded (distance < MIN_DISTANCE)
    let mut mol = Molecule::new("overlap");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
    infer_bonds(&mut mol, None).unwrap();
    assert_eq!(mol.bond_count(), 0);
}
