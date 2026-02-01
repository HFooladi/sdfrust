//! Coordinate transformations for molecular structures.

use crate::Molecule;
use nalgebra::{Matrix3, Rotation3, Unit, Vector3};

/// Rotate a molecule around an axis by a given angle.
///
/// Uses Rodrigues' rotation formula via nalgebra. The rotation is performed
/// around the origin - center the molecule first if rotation around the
/// centroid is desired.
///
/// # Arguments
///
/// * `molecule` - The molecule to rotate (modified in place)
/// * `axis` - The rotation axis as [x, y, z] (will be normalized)
/// * `angle` - The rotation angle in radians
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::rotate;
/// use std::f64::consts::PI;
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
///
/// // Rotate 90 degrees around Z axis
/// rotate(&mut mol, [0.0, 0.0, 1.0], PI / 2.0);
///
/// // Point (1, 0, 0) rotated 90° around Z should be near (0, 1, 0)
/// assert!((mol.atoms[0].x - 0.0).abs() < 1e-10);
/// assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
/// ```
pub fn rotate(molecule: &mut Molecule, axis: [f64; 3], angle: f64) {
    if molecule.atoms.is_empty() {
        return;
    }

    let axis_vec = Vector3::new(axis[0], axis[1], axis[2]);
    if axis_vec.norm() < 1e-10 {
        return; // Zero axis, no rotation
    }

    let unit_axis = Unit::new_normalize(axis_vec);
    let rotation = Rotation3::from_axis_angle(&unit_axis, angle);

    for atom in &mut molecule.atoms {
        let point = Vector3::new(atom.x, atom.y, atom.z);
        let rotated = rotation * point;
        atom.x = rotated.x;
        atom.y = rotated.y;
        atom.z = rotated.z;
    }
}

/// Apply a 3x3 rotation matrix to a molecule.
///
/// The matrix should be a proper rotation matrix (orthogonal with determinant 1).
/// No validation is performed - using a non-rotation matrix will give
/// unexpected results.
///
/// # Arguments
///
/// * `molecule` - The molecule to transform (modified in place)
/// * `matrix` - A 3x3 rotation matrix in row-major order
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::apply_rotation_matrix;
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
///
/// // Identity matrix - no change
/// let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// apply_rotation_matrix(&mut mol, &identity);
///
/// assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
/// ```
pub fn apply_rotation_matrix(molecule: &mut Molecule, matrix: &[[f64; 3]; 3]) {
    if molecule.atoms.is_empty() {
        return;
    }

    let rot = Matrix3::new(
        matrix[0][0],
        matrix[0][1],
        matrix[0][2],
        matrix[1][0],
        matrix[1][1],
        matrix[1][2],
        matrix[2][0],
        matrix[2][1],
        matrix[2][2],
    );

    for atom in &mut molecule.atoms {
        let point = Vector3::new(atom.x, atom.y, atom.z);
        let transformed = rot * point;
        atom.x = transformed.x;
        atom.y = transformed.y;
        atom.z = transformed.z;
    }
}

/// Apply a rotation matrix and translation to a molecule.
///
/// First applies the rotation, then the translation.
/// This is equivalent to the affine transformation: p' = R*p + t
///
/// # Arguments
///
/// * `molecule` - The molecule to transform (modified in place)
/// * `rotation` - A 3x3 rotation matrix in row-major order
/// * `translation` - A translation vector [dx, dy, dz]
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::apply_transform;
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
///
/// // Identity rotation + translation by (1, 2, 3)
/// let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
/// apply_transform(&mut mol, &identity, [1.0, 2.0, 3.0]);
///
/// assert!((mol.atoms[0].x - 2.0).abs() < 1e-10);
/// assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
/// assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
/// ```
pub fn apply_transform(molecule: &mut Molecule, rotation: &[[f64; 3]; 3], translation: [f64; 3]) {
    if molecule.atoms.is_empty() {
        return;
    }

    let rot = Matrix3::new(
        rotation[0][0],
        rotation[0][1],
        rotation[0][2],
        rotation[1][0],
        rotation[1][1],
        rotation[1][2],
        rotation[2][0],
        rotation[2][1],
        rotation[2][2],
    );

    for atom in &mut molecule.atoms {
        let point = Vector3::new(atom.x, atom.y, atom.z);
        let transformed = rot * point;
        atom.x = transformed.x + translation[0];
        atom.y = transformed.y + translation[1];
        atom.z = transformed.z + translation[2];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;
    use std::f64::consts::PI;

    #[test]
    fn test_rotate_identity() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

        // Zero angle rotation
        rotate(&mut mol, [0.0, 0.0, 1.0], 0.0);

        assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_90_z_axis() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

        rotate(&mut mol, [0.0, 0.0, 1.0], PI / 2.0);

        // (1, 0, 0) rotated 90° around Z -> (0, 1, 0)
        assert!((mol.atoms[0].x).abs() < 1e-10);
        assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].z).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_180_x_axis() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 1.0, 0.0));

        rotate(&mut mol, [1.0, 0.0, 0.0], PI);

        // (0, 1, 0) rotated 180° around X -> (0, -1, 0)
        assert!((mol.atoms[0].x).abs() < 1e-10);
        assert!((mol.atoms[0].y + 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].z).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_preserves_distances() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));

        let d01_before = mol.atoms[0].distance_to(&mol.atoms[1]);
        let d02_before = mol.atoms[0].distance_to(&mol.atoms[2]);
        let d12_before = mol.atoms[1].distance_to(&mol.atoms[2]);

        // Arbitrary rotation
        rotate(&mut mol, [1.0, 2.0, 3.0], 1.234);

        let d01_after = mol.atoms[0].distance_to(&mol.atoms[1]);
        let d02_after = mol.atoms[0].distance_to(&mol.atoms[2]);
        let d12_after = mol.atoms[1].distance_to(&mol.atoms[2]);

        assert!((d01_before - d01_after).abs() < 1e-10);
        assert!((d02_before - d02_after).abs() < 1e-10);
        assert!((d12_before - d12_after).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_empty_molecule() {
        let mut mol = Molecule::new("empty");
        rotate(&mut mol, [0.0, 0.0, 1.0], PI / 2.0);
        assert!(mol.atoms.is_empty());
    }

    #[test]
    fn test_rotate_zero_axis() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

        // Zero axis should not change anything
        rotate(&mut mol, [0.0, 0.0, 0.0], PI / 2.0);

        assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_rotation_matrix_identity() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        apply_rotation_matrix(&mut mol, &identity);

        assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_transform() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        apply_transform(&mut mol, &identity, [5.0, 10.0, 15.0]);

        assert!((mol.atoms[0].x - 6.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 10.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_apply_transform_with_rotation() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

        // 90° rotation around Z, then translate by (0, 0, 5)
        let rot_90_z = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        apply_transform(&mut mol, &rot_90_z, [0.0, 0.0, 5.0]);

        // (1, 0, 0) -> (0, 1, 0) -> (0, 1, 5)
        assert!((mol.atoms[0].x).abs() < 1e-10);
        assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 5.0).abs() < 1e-10);
    }
}
