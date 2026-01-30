//! Tests for molecular descriptors.
//!
//! These tests verify the descriptor calculations against known values
//! from PubChem and calculated reference data.

use sdfrust::{
    parse_sdf_file, Atom, Bond, BondOrder, Molecule,
    descriptors::{
        atomic_weight, bond_type_counts, exact_mass, get_element, heavy_atom_count,
        molecular_weight, monoisotopic_mass, ring_atoms, ring_bonds, ring_count,
        rotatable_bond_count,
    },
};

// ============================================================
// Element Data Tests
// ============================================================

#[test]
fn test_element_data_carbon() {
    let c = get_element("C").unwrap();
    assert_eq!(c.atomic_number, 6);
    assert_eq!(c.symbol, "C");
    // IUPAC 2021 standard atomic weight
    assert!((c.atomic_weight - 12.011).abs() < 0.001);
    // Carbon-12 is exactly 12.0 by definition
    assert!((c.monoisotopic_mass - 12.0).abs() < 0.0001);
}

#[test]
fn test_element_data_common_elements() {
    // Verify common organic elements
    let elements = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"];
    for sym in elements {
        assert!(get_element(sym).is_some(), "Missing element: {}", sym);
    }
}

#[test]
fn test_element_data_metals() {
    // Verify common metals
    let metals = ["Na", "K", "Ca", "Mg", "Fe", "Cu", "Zn"];
    for sym in metals {
        assert!(get_element(sym).is_some(), "Missing metal: {}", sym);
    }
}

#[test]
fn test_atomic_weight_helper() {
    assert!((atomic_weight("O").unwrap() - 15.999).abs() < 0.001);
    assert!((atomic_weight("N").unwrap() - 14.007).abs() < 0.001);
    assert!(atomic_weight("Xx").is_none());
}

#[test]
fn test_monoisotopic_mass_helper() {
    // Carbon-12 is exactly 12.0
    assert!((monoisotopic_mass("C").unwrap() - 12.0).abs() < 0.0001);
    // Oxygen-16
    assert!((monoisotopic_mass("O").unwrap() - 15.99491).abs() < 0.0001);
}

#[test]
fn test_element_case_insensitive() {
    assert!(get_element("c").is_some());
    assert!(get_element("C").is_some());
    assert!(get_element("cl").is_some());
    assert!(get_element("Cl").is_some());
    assert!(get_element("CL").is_some());
}

// ============================================================
// Molecular Weight Tests
// ============================================================

#[test]
fn test_molecular_weight_water() {
    let mut mol = Molecule::new("water");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));

    let mw = molecular_weight(&mol).unwrap();
    // H2O: 2*1.008 + 15.999 = 18.015
    assert!((mw - 18.015).abs() < 0.001);
}

#[test]
fn test_molecular_weight_methane() {
    let mut mol = Molecule::new("methane");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.63, 0.63, 0.63));
    mol.atoms.push(Atom::new(2, "H", -0.63, -0.63, 0.63));
    mol.atoms.push(Atom::new(3, "H", -0.63, 0.63, -0.63));
    mol.atoms.push(Atom::new(4, "H", 0.63, -0.63, -0.63));

    let mw = molecular_weight(&mol).unwrap();
    // CH4: 12.011 + 4*1.008 = 16.043
    assert!((mw - 16.043).abs() < 0.001);
}

#[test]
fn test_molecular_weight_aspirin() {
    // Load aspirin from test file
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let mw = molecular_weight(&mol).unwrap();
    // PubChem value: 180.16
    assert!((mw - 180.16).abs() < 0.05, "Aspirin MW: expected ~180.16, got {}", mw);
}

#[test]
fn test_molecular_weight_caffeine() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    let mw = molecular_weight(&mol).unwrap();
    // PubChem value: 194.19
    assert!((mw - 194.19).abs() < 0.05, "Caffeine MW: expected ~194.19, got {}", mw);
}

#[test]
fn test_molecular_weight_unknown_element() {
    let mut mol = Molecule::new("unknown");
    mol.atoms.push(Atom::new(0, "Xx", 0.0, 0.0, 0.0));
    assert!(molecular_weight(&mol).is_none());
}

#[test]
fn test_molecular_weight_empty() {
    let mol = Molecule::new("empty");
    let mw = molecular_weight(&mol).unwrap();
    assert!((mw - 0.0).abs() < 0.001);
}

#[test]
fn test_molecular_weight_via_molecule_method() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let mw = mol.molecular_weight().unwrap();
    assert!((mw - 180.16).abs() < 0.05);
}

// ============================================================
// Exact Mass Tests
// ============================================================

#[test]
fn test_exact_mass_water() {
    let mut mol = Molecule::new("water");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));

    let mass = exact_mass(&mol).unwrap();
    // H2O: 2*1.00783 + 15.99491 = 18.01056
    assert!((mass - 18.01056).abs() < 0.001);
}

#[test]
fn test_exact_mass_aspirin() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let mass = exact_mass(&mol).unwrap();
    // PubChem value: 180.04225873
    assert!((mass - 180.042).abs() < 0.01, "Aspirin exact mass: expected ~180.042, got {}", mass);
}

#[test]
fn test_exact_mass_caffeine() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    let mass = exact_mass(&mol).unwrap();
    // PubChem value: 194.08037557
    assert!((mass - 194.080).abs() < 0.01, "Caffeine exact mass: expected ~194.080, got {}", mass);
}

#[test]
fn test_exact_mass_via_molecule_method() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let mass = mol.exact_mass().unwrap();
    assert!((mass - 180.042).abs() < 0.01);
}

// ============================================================
// Heavy Atom Count Tests
// ============================================================

#[test]
fn test_heavy_atom_count_water() {
    let mut mol = Molecule::new("water");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));

    assert_eq!(heavy_atom_count(&mol), 1);
}

#[test]
fn test_heavy_atom_count_aspirin() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    // PubChem value: 13
    assert_eq!(heavy_atom_count(&mol), 13);
}

#[test]
fn test_heavy_atom_count_caffeine() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    // PubChem value: 14
    assert_eq!(heavy_atom_count(&mol), 14);
}

#[test]
fn test_heavy_atom_count_empty() {
    let mol = Molecule::new("empty");
    assert_eq!(heavy_atom_count(&mol), 0);
}

#[test]
fn test_heavy_atom_count_all_hydrogen() {
    let mut mol = Molecule::new("h2");
    mol.atoms.push(Atom::new(0, "H", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.74, 0.0, 0.0));
    assert_eq!(heavy_atom_count(&mol), 0);
}

#[test]
fn test_heavy_atom_count_via_molecule_method() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    assert_eq!(mol.heavy_atom_count(), 13);
}

// ============================================================
// Bond Type Count Tests
// ============================================================

#[test]
fn test_bond_type_counts_water() {
    let mut mol = Molecule::new("water");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

    let counts = bond_type_counts(&mol);
    assert_eq!(counts.get(&BondOrder::Single), Some(&2));
    assert_eq!(counts.get(&BondOrder::Double), None);
}

#[test]
fn test_bond_type_counts_co2() {
    let mut mol = Molecule::new("co2");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "O", -1.16, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "O", 1.16, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
    mol.bonds.push(Bond::new(0, 2, BondOrder::Double));

    let counts = bond_type_counts(&mol);
    assert_eq!(counts.get(&BondOrder::Single), None);
    assert_eq!(counts.get(&BondOrder::Double), Some(&2));
}

#[test]
fn test_bond_type_counts_benzene_aromatic() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    let counts = bond_type_counts(&mol);
    assert_eq!(counts.get(&BondOrder::Aromatic), Some(&6));
}

#[test]
fn test_bond_type_counts_empty() {
    let mol = Molecule::new("empty");
    let counts = bond_type_counts(&mol);
    assert!(counts.is_empty());
}

#[test]
fn test_bond_type_counts_via_molecule_method() {
    let mut mol = Molecule::new("ethene");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

    let counts = mol.bond_type_counts();
    assert_eq!(counts.get(&BondOrder::Double), Some(&1));
}

// ============================================================
// Ring Count Tests
// ============================================================

#[test]
fn test_ring_count_benzene() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        let angle = std::f64::consts::PI * 2.0 * i as f64 / 6.0;
        mol.atoms.push(Atom::new(i, "C", angle.cos(), angle.sin(), 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    assert_eq!(ring_count(&mol), 1);
}

#[test]
fn test_ring_count_naphthalene() {
    // Two fused benzene rings = 2 independent cycles
    let mut mol = Molecule::new("naphthalene");
    for i in 0..10 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    // First ring: 0-1-2-3-4-5-0
    mol.bonds.push(Bond::new(0, 1, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(1, 2, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(2, 3, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(3, 4, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(4, 5, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(5, 0, BondOrder::Aromatic));
    // Second ring shares edge 3-4: 3-4-6-7-8-9-3
    mol.bonds.push(Bond::new(4, 6, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(6, 7, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(7, 8, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(8, 9, BondOrder::Aromatic));
    mol.bonds.push(Bond::new(9, 3, BondOrder::Aromatic));

    assert_eq!(ring_count(&mol), 2);
}

#[test]
fn test_ring_count_propane() {
    let mut mol = Molecule::new("propane");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "C", 3.0, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 2, BondOrder::Single));

    assert_eq!(ring_count(&mol), 0);
}

#[test]
fn test_ring_count_empty() {
    let mol = Molecule::new("empty");
    assert_eq!(ring_count(&mol), 0);
}

#[test]
fn test_ring_count_cyclopropane() {
    let mut mol = Molecule::new("cyclopropane");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "C", 0.5, 0.87, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
    mol.bonds.push(Bond::new(2, 0, BondOrder::Single));

    assert_eq!(ring_count(&mol), 1);
}

#[test]
fn test_ring_count_via_molecule_method() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    assert_eq!(mol.ring_count(), 1);
}

#[test]
fn test_ring_count_caffeine() {
    // Caffeine has 2 fused rings (purine scaffold)
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    assert_eq!(ring_count(&mol), 2, "Caffeine should have 2 rings (purine scaffold)");
}

// ============================================================
// Ring Membership Tests
// ============================================================

#[test]
fn test_ring_atoms_benzene() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    let in_ring = ring_atoms(&mol);
    assert_eq!(in_ring.len(), 6);
    assert!(in_ring.iter().all(|&r| r), "All benzene atoms should be in ring");
}

#[test]
fn test_ring_atoms_toluene_like() {
    // Ring with one substituent
    let mut mol = Molecule::new("toluene_like");
    // Benzene ring: atoms 0-5
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    // Methyl: atom 6
    mol.atoms.push(Atom::new(6, "C", 2.0, 0.0, 0.0));

    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }
    mol.bonds.push(Bond::new(0, 6, BondOrder::Single));

    let in_ring = ring_atoms(&mol);
    assert_eq!(in_ring.len(), 7);
    // Ring carbons (0-5) should be in ring
    for i in 0..6 {
        assert!(in_ring[i], "Atom {} should be in ring", i);
    }
    // Methyl (6) should not be in ring
    assert!(!in_ring[6], "Methyl carbon should not be in ring");
}

#[test]
fn test_ring_bonds_benzene() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    let in_ring = ring_bonds(&mol);
    assert_eq!(in_ring.len(), 6);
    assert!(in_ring.iter().all(|&r| r), "All benzene bonds should be in ring");
}

#[test]
fn test_ring_bonds_toluene_like() {
    let mut mol = Molecule::new("toluene_like");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    mol.atoms.push(Atom::new(6, "C", 2.0, 0.0, 0.0));

    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }
    mol.bonds.push(Bond::new(0, 6, BondOrder::Single)); // Bond 6: exocyclic

    let in_ring = ring_bonds(&mol);
    assert_eq!(in_ring.len(), 7);
    // Ring bonds (0-5) should be in ring
    for i in 0..6 {
        assert!(in_ring[i], "Bond {} should be in ring", i);
    }
    // Exocyclic bond (6) should not be in ring
    assert!(!in_ring[6], "Exocyclic bond should not be in ring");
}

#[test]
fn test_is_atom_in_ring_via_molecule_method() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    for i in 0..6 {
        assert!(mol.is_atom_in_ring(i));
    }
    // Out of bounds should return false
    assert!(!mol.is_atom_in_ring(100));
}

#[test]
fn test_is_bond_in_ring_via_molecule_method() {
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    for i in 0..6 {
        assert!(mol.is_bond_in_ring(i));
    }
    // Out of bounds should return false
    assert!(!mol.is_bond_in_ring(100));
}

// ============================================================
// Rotatable Bond Tests
// ============================================================

#[test]
fn test_rotatable_bond_count_benzene() {
    // Benzene has no rotatable bonds (all aromatic, in ring)
    let mut mol = Molecule::new("benzene");
    for i in 0..6 {
        mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
    }
    for i in 0..6 {
        mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
    }

    assert_eq!(rotatable_bond_count(&mol), 0);
}

#[test]
fn test_rotatable_bond_count_ethane() {
    // Ethane: C-C with hydrogens (CH3-CH3)
    // Both carbons have heavy_degree = 1 (only bonded to one heavy atom - each other)
    // Therefore the C-C bond is NOT rotatable (both ends are terminal)
    let mut mol = Molecule::new("ethane");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    // 6 hydrogens
    for i in 2..8 {
        mol.atoms.push(Atom::new(i, "H", 0.0, 0.0, 0.0));
    }
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single)); // C-C
    mol.bonds.push(Bond::new(0, 2, BondOrder::Single)); // C-H
    mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 5, BondOrder::Single)); // C-H
    mol.bonds.push(Bond::new(1, 6, BondOrder::Single));
    mol.bonds.push(Bond::new(1, 7, BondOrder::Single));

    // C-C bond is NOT rotatable because both carbons are terminal
    assert_eq!(rotatable_bond_count(&mol), 0);
}

#[test]
fn test_rotatable_bond_count_aspirin() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    // PubChem value: 3
    // RDKit definition: single, non-ring, both atoms non-terminal (heavy_degree > 1)
    assert_eq!(rotatable_bond_count(&mol), 3, "Aspirin should have 3 rotatable bonds");
}

#[test]
fn test_rotatable_bond_count_caffeine() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    // PubChem value: 0 (rigid fused ring system)
    assert_eq!(rotatable_bond_count(&mol), 0, "Caffeine should have 0 rotatable bonds");
}

#[test]
fn test_rotatable_bond_count_empty() {
    let mol = Molecule::new("empty");
    assert_eq!(rotatable_bond_count(&mol), 0);
}

#[test]
fn test_rotatable_bond_count_single_bond_terminal() {
    // Two atoms connected by single bond - both terminal
    let mut mol = Molecule::new("simple");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

    // Both atoms have degree 1 (terminal), so not rotatable
    assert_eq!(rotatable_bond_count(&mol), 0);
}

#[test]
fn test_rotatable_bond_count_via_molecule_method() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    assert_eq!(mol.rotatable_bond_count(), 3);
}

// ============================================================
// Edge Case Tests
// ============================================================

#[test]
fn test_descriptors_single_atom() {
    let mut mol = Molecule::new("carbon");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

    assert!((molecular_weight(&mol).unwrap() - 12.011).abs() < 0.001);
    assert!((exact_mass(&mol).unwrap() - 12.0).abs() < 0.001);
    assert_eq!(heavy_atom_count(&mol), 1);
    assert_eq!(ring_count(&mol), 0);
    assert_eq!(rotatable_bond_count(&mol), 0);
}

#[test]
fn test_descriptors_disconnected_fragments() {
    // Two separate fragments
    let mut mol = Molecule::new("two_fragments");
    // Fragment 1: C-C
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    // Fragment 2: O (isolated)
    mol.atoms.push(Atom::new(2, "O", 5.0, 0.0, 0.0));

    // MW should include all atoms
    let mw = molecular_weight(&mol).unwrap();
    // 2*12.011 + 15.999 = 40.021
    assert!((mw - 40.021).abs() < 0.001);

    // Heavy atom count should be 3
    assert_eq!(heavy_atom_count(&mol), 3);

    // No rings
    assert_eq!(ring_count(&mol), 0);
}

#[test]
fn test_descriptors_with_deuterium() {
    let mut mol = Molecule::new("heavy_water");
    mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "D", 0.96, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "D", -0.24, 0.93, 0.0));

    let mw = molecular_weight(&mol).unwrap();
    // D2O: 2*2.014 + 15.999 ≈ 20.027
    assert!((mw - 20.027).abs() < 0.01, "D2O MW: expected ~20.027, got {}", mw);

    // Deuterium counts as hydrogen for heavy atom count
    assert_eq!(heavy_atom_count(&mol), 1);
}

// ============================================================
// Integration Tests with Real Files
// ============================================================

#[test]
fn test_glucose_descriptors() {
    let mol = parse_sdf_file("tests/test_data/glucose.sdf").unwrap();

    // Glucose (actually sucrose in the file based on name) - check formula
    let formula = mol.formula();
    println!("Glucose file formula: {}", formula);

    // Should have reasonable MW
    let mw = mol.molecular_weight().unwrap();
    assert!(mw > 100.0, "Glucose MW should be > 100");

    // Should have heavy atoms
    let heavy = mol.heavy_atom_count();
    assert!(heavy > 0, "Should have heavy atoms");
}

#[test]
fn test_methionine_descriptors() {
    let mol = parse_sdf_file("tests/test_data/methionine.sdf").unwrap();

    // Methionine: C5H11NO2S
    // MW ≈ 149.21
    let mw = mol.molecular_weight().unwrap();
    assert!((mw - 149.2).abs() < 1.0, "Methionine MW: expected ~149.2, got {}", mw);

    // Contains S
    let has_sulfur = mol.atoms.iter().any(|a| a.element == "S");
    assert!(has_sulfur, "Methionine should contain sulfur");
}

#[test]
fn test_galactose_descriptors() {
    let mol = parse_sdf_file("tests/test_data/galactose.sdf").unwrap();

    // Galactose: C6H12O6
    // MW ≈ 180.16
    let mw = mol.molecular_weight().unwrap();
    assert!((mw - 180.16).abs() < 1.0, "Galactose MW: expected ~180.16, got {}", mw);
}
