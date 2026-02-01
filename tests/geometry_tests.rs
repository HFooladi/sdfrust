//! Integration tests for the geometry module.
//!
//! These tests verify the geometry operations work correctly with
//! the overall library structure.

#![cfg(feature = "geometry")]

use std::f64::consts::PI;

use sdfrust::{Atom, Bond, BondOrder, Molecule};

// ============================================================
// Distance Matrix Tests
// ============================================================

#[test]
fn test_distance_matrix_empty_molecule() {
    let mol = Molecule::new("empty");
    let matrix = mol.distance_matrix();
    assert!(matrix.is_empty());
}

#[test]
fn test_distance_matrix_single_atom() {
    let mut mol = Molecule::new("single");
    mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

    let matrix = mol.distance_matrix();
    assert_eq!(matrix.len(), 1);
    assert_eq!(matrix[0].len(), 1);
    assert!((matrix[0][0]).abs() < 1e-10);
}

#[test]
fn test_distance_matrix_two_atoms() {
    let mut mol = Molecule::new("two");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 3.0, 4.0, 0.0));

    let matrix = mol.distance_matrix();
    assert_eq!(matrix.len(), 2);
    assert!((matrix[0][1] - 5.0).abs() < 1e-10);
    assert!((matrix[1][0] - 5.0).abs() < 1e-10);
}

#[test]
fn test_distance_matrix_symmetric() {
    let mut mol = Molecule::new("triangle");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "C", 0.5, 0.866, 0.0));

    let matrix = mol.distance_matrix();

    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (matrix[i][j] - matrix[j][i]).abs() < 1e-10,
                "Matrix should be symmetric"
            );
        }
    }
}

#[test]
fn test_distance_matrix_correct_values() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
    mol.atoms.push(Atom::new(3, "C", 0.0, 0.0, 1.0));

    let matrix = mol.distance_matrix();

    // Diagonal should be zero
    for i in 0..4 {
        assert!((matrix[i][i]).abs() < 1e-10);
    }

    // Distance 0-1, 0-2, 0-3 should all be 1.0
    assert!((matrix[0][1] - 1.0).abs() < 1e-10);
    assert!((matrix[0][2] - 1.0).abs() < 1e-10);
    assert!((matrix[0][3] - 1.0).abs() < 1e-10);

    // Distance 1-2 should be sqrt(2)
    assert!((matrix[1][2] - 2.0_f64.sqrt()).abs() < 1e-10);
}

// ============================================================
// RMSD Tests
// ============================================================

#[test]
fn test_rmsd_identical_molecules() {
    let mut mol1 = Molecule::new("mol1");
    mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol1.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    mol1.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));

    let mol2 = mol1.clone();

    let rmsd = mol1.rmsd_to(&mol2).unwrap();
    assert!(
        rmsd.abs() < 1e-10,
        "RMSD of identical molecules should be 0"
    );
}

#[test]
fn test_rmsd_translated_molecule() {
    let mut mol1 = Molecule::new("mol1");
    mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol1.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));

    let mut mol2 = mol1.clone();
    mol2.translate(1.0, 0.0, 0.0);

    let rmsd = mol1.rmsd_to(&mol2).unwrap();
    assert!(
        (rmsd - 1.0).abs() < 1e-10,
        "RMSD should be 1.0 for 1A translation"
    );
}

#[test]
fn test_rmsd_atom_count_mismatch() {
    let mut mol1 = Molecule::new("mol1");
    mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

    let mut mol2 = Molecule::new("mol2");
    mol2.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol2.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));

    let result = mol1.rmsd_to(&mol2);
    assert!(
        result.is_err(),
        "RMSD should fail for different atom counts"
    );
}

#[test]
fn test_rmsd_empty_molecules() {
    let mol1 = Molecule::new("mol1");
    let mol2 = Molecule::new("mol2");

    let rmsd = mol1.rmsd_to(&mol2).unwrap();
    assert!(rmsd.abs() < 1e-10, "RMSD of empty molecules should be 0");
}

// ============================================================
// Rotation Tests
// ============================================================

#[test]
fn test_rotate_identity() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

    // Zero angle rotation
    mol.rotate([0.0, 0.0, 1.0], 0.0);

    assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
    assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
    assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
}

#[test]
fn test_rotate_90_degrees_z_axis() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

    mol.rotate([0.0, 0.0, 1.0], PI / 2.0);

    // (1, 0, 0) rotated 90° around Z -> (0, 1, 0)
    assert!(mol.atoms[0].x.abs() < 1e-10, "x should be 0");
    assert!((mol.atoms[0].y - 1.0).abs() < 1e-10, "y should be 1");
    assert!(mol.atoms[0].z.abs() < 1e-10, "z should be 0");
}

#[test]
fn test_rotate_180_degrees_x_axis() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 1.0, 0.0));

    mol.rotate([1.0, 0.0, 0.0], PI);

    // (0, 1, 0) rotated 180° around X -> (0, -1, 0)
    assert!(mol.atoms[0].x.abs() < 1e-10);
    assert!((mol.atoms[0].y + 1.0).abs() < 1e-10);
    assert!(mol.atoms[0].z.abs() < 1e-10);
}

#[test]
fn test_rotate_preserves_distances() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
    mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

    let d01_before = mol.atoms[0].distance_to(&mol.atoms[1]);
    let d02_before = mol.atoms[0].distance_to(&mol.atoms[2]);
    let d12_before = mol.atoms[1].distance_to(&mol.atoms[2]);

    // Arbitrary rotation
    mol.rotate([1.0, 2.0, 3.0], 1.234);

    let d01_after = mol.atoms[0].distance_to(&mol.atoms[1]);
    let d02_after = mol.atoms[0].distance_to(&mol.atoms[2]);
    let d12_after = mol.atoms[1].distance_to(&mol.atoms[2]);

    assert!(
        (d01_before - d01_after).abs() < 1e-10,
        "Rotation should preserve distances"
    );
    assert!((d02_before - d02_after).abs() < 1e-10);
    assert!((d12_before - d12_after).abs() < 1e-10);
}

#[test]
fn test_rotate_empty_molecule() {
    let mut mol = Molecule::new("empty");
    mol.rotate([0.0, 0.0, 1.0], PI / 2.0);
    assert!(mol.atoms.is_empty());
}

#[test]
fn test_rotate_arbitrary_axis() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "C", 0.0, 1.0, 0.0));

    let d_before = mol.atoms[0].distance_to(&mol.atoms[1]);

    // Rotate around (1, 1, 1) axis
    mol.rotate([1.0, 1.0, 1.0], 2.0 * PI / 3.0);

    let d_after = mol.atoms[0].distance_to(&mol.atoms[1]);
    assert!(
        (d_before - d_after).abs() < 1e-10,
        "Distance should be preserved"
    );
}

// ============================================================
// Transform Tests
// ============================================================

#[test]
fn test_apply_rotation_matrix_identity() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    mol.apply_rotation_matrix(&identity);

    assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
    assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
    assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
}

#[test]
fn test_apply_transform_rotation_and_translation() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

    // 90° rotation around Z, then translate by (0, 0, 5)
    let rot_90_z = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    mol.apply_transform(&rot_90_z, [0.0, 0.0, 5.0]);

    // (1, 0, 0) -> (0, 1, 0) -> (0, 1, 5)
    assert!(mol.atoms[0].x.abs() < 1e-10);
    assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
    assert!((mol.atoms[0].z - 5.0).abs() < 1e-10);
}

// ============================================================
// Integration Tests with Real Molecules
// ============================================================

#[test]
fn test_geometry_with_methane() {
    let mut mol = Molecule::new("methane");

    // Tetrahedral methane
    mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    mol.atoms.push(Atom::new(1, "H", 1.089, 0.0, 0.0));
    mol.atoms.push(Atom::new(2, "H", -0.363, 1.027, 0.0));
    mol.atoms.push(Atom::new(3, "H", -0.363, -0.513, 0.889));
    mol.atoms.push(Atom::new(4, "H", -0.363, -0.513, -0.889));

    // Test distance matrix
    let matrix = mol.distance_matrix();
    assert_eq!(matrix.len(), 5);

    // All C-H distances should be similar (~1.089)
    for i in 1..5 {
        assert!(
            (matrix[0][i] - 1.089).abs() < 0.01,
            "C-H distance should be ~1.089"
        );
    }

    // Test RMSD after rotation
    let mol_original = mol.clone();
    mol.rotate([1.0, 1.0, 1.0], PI / 4.0);

    // RMSD should be non-zero after rotation (not aligned)
    let rmsd = mol.rmsd_to(&mol_original).unwrap();
    assert!(rmsd > 0.0, "RMSD should be non-zero after rotation");

    // But distances should be preserved
    let matrix_after = mol.distance_matrix();
    for i in 0..5 {
        for j in 0..5 {
            assert!(
                (matrix[i][j] - matrix_after[i][j]).abs() < 1e-10,
                "Distances should be preserved after rotation"
            );
        }
    }
}

#[test]
fn test_center_and_rotate() {
    let mut mol = Molecule::new("test");
    mol.atoms.push(Atom::new(0, "C", 10.0, 10.0, 10.0));
    mol.atoms.push(Atom::new(1, "C", 11.0, 10.0, 10.0));

    // Center first
    mol.center();
    let centroid = mol.centroid().unwrap();
    assert!(centroid.0.abs() < 1e-10);
    assert!(centroid.1.abs() < 1e-10);
    assert!(centroid.2.abs() < 1e-10);

    // Then rotate
    mol.rotate([0.0, 0.0, 1.0], PI / 2.0);

    // Centroid should still be at origin
    let centroid_after = mol.centroid().unwrap();
    assert!(
        centroid_after.0.abs() < 1e-10,
        "Centroid x should be 0 after rotation"
    );
    assert!(centroid_after.1.abs() < 1e-10);
    assert!(centroid_after.2.abs() < 1e-10);
}
