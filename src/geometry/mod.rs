//! Geometric operations for molecular structures.
//!
//! This module provides optional nalgebra-based geometric operations including
//! rotation, distance matrix calculation, and RMSD computation.
//!
//! Enable with: `sdfrust = { features = ["geometry"] }`
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom};
//! use sdfrust::geometry;
//! use std::f64::consts::PI;
//!
//! let mut mol = Molecule::new("example");
//! mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 2.0, 0.0, 0.0));
//!
//! // Get distance matrix
//! let matrix = mol.distance_matrix();
//! assert!((matrix[0][1] - 1.0).abs() < 1e-10);
//!
//! // Rotate around Z axis
//! mol.rotate([0.0, 0.0, 1.0], PI / 2.0);
//! ```

mod distance;
mod rmsd;
mod transform;

pub use distance::distance_matrix;
pub use rmsd::{rmsd_from_coords, rmsd_to};
pub use transform::{apply_rotation_matrix, apply_transform, rotate};

use crate::{Molecule, SdfError};

// Extension methods for Molecule
impl Molecule {
    /// Rotate the molecule around an axis by a given angle.
    ///
    /// Uses Rodrigues' rotation formula. The rotation is performed around
    /// the origin - center the molecule first if rotation around the
    /// centroid is desired.
    ///
    /// # Arguments
    ///
    /// * `axis` - The rotation axis as [x, y, z] (will be normalized)
    /// * `angle` - The rotation angle in radians
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    /// use std::f64::consts::PI;
    ///
    /// let mut mol = Molecule::new("test");
    /// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
    ///
    /// // Rotate 90 degrees around Z axis
    /// mol.rotate([0.0, 0.0, 1.0], PI / 2.0);
    ///
    /// // Point (1, 0, 0) rotated 90° around Z should be near (0, 1, 0)
    /// assert!((mol.atoms[0].x - 0.0).abs() < 1e-10);
    /// assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
    /// ```
    pub fn rotate(&mut self, axis: [f64; 3], angle: f64) {
        rotate(self, axis, angle);
    }

    /// Apply a 3x3 rotation matrix to the molecule.
    ///
    /// The matrix should be a proper rotation matrix (orthogonal with determinant 1).
    /// No validation is performed.
    ///
    /// # Arguments
    ///
    /// * `matrix` - A 3x3 rotation matrix in row-major order
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    ///
    /// let mut mol = Molecule::new("test");
    /// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
    ///
    /// // Identity matrix
    /// let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// mol.apply_rotation_matrix(&identity);
    /// ```
    pub fn apply_rotation_matrix(&mut self, matrix: &[[f64; 3]; 3]) {
        apply_rotation_matrix(self, matrix);
    }

    /// Apply a rotation matrix and translation to the molecule.
    ///
    /// First applies the rotation, then the translation.
    /// This is equivalent to the affine transformation: p' = R*p + t
    ///
    /// # Arguments
    ///
    /// * `rotation` - A 3x3 rotation matrix in row-major order
    /// * `translation` - A translation vector [dx, dy, dz]
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    ///
    /// let mut mol = Molecule::new("test");
    /// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
    ///
    /// let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// mol.apply_transform(&identity, [1.0, 2.0, 3.0]);
    /// ```
    pub fn apply_transform(&mut self, rotation: &[[f64; 3]; 3], translation: [f64; 3]) {
        apply_transform(self, rotation, translation);
    }

    /// Compute the pairwise distance matrix for all atoms.
    ///
    /// Returns an NxN matrix where entry [i][j] is the Euclidean distance
    /// between atom i and atom j in Angstroms.
    ///
    /// # Returns
    ///
    /// A symmetric matrix of pairwise distances. The diagonal is all zeros.
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    ///
    /// let mut mol = Molecule::new("test");
    /// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    /// mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
    ///
    /// let matrix = mol.distance_matrix();
    /// assert!((matrix[0][1] - 1.0).abs() < 1e-10);
    /// ```
    pub fn distance_matrix(&self) -> Vec<Vec<f64>> {
        distance_matrix(self)
    }

    /// Calculate RMSD to another molecule.
    ///
    /// Computes the root mean square deviation of atomic positions.
    /// The molecules must have the same number of atoms.
    /// No alignment is performed - atoms are compared directly by index.
    ///
    /// # Arguments
    ///
    /// * `other` - The other molecule to compare to
    ///
    /// # Returns
    ///
    /// The RMSD value in Angstroms, or an error if the molecules have different atom counts.
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    ///
    /// let mut mol1 = Molecule::new("mol1");
    /// mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
    ///
    /// let mut mol2 = Molecule::new("mol2");
    /// mol2.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
    ///
    /// let rmsd = mol1.rmsd_to(&mol2).unwrap();
    /// assert!((rmsd - 1.0).abs() < 1e-10);
    /// ```
    pub fn rmsd_to(&self, other: &Molecule) -> Result<f64, SdfError> {
        rmsd_to(self, other)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;
    use std::f64::consts::PI;

    #[test]
    fn test_molecule_rotate_method() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));

        mol.rotate([0.0, 0.0, 1.0], PI / 2.0);

        assert!((mol.atoms[0].x).abs() < 1e-10);
        assert!((mol.atoms[0].y - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_molecule_distance_matrix_method() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 3.0, 4.0, 0.0));

        let matrix = mol.distance_matrix();
        assert!((matrix[0][1] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_molecule_rmsd_to_method() {
        let mut mol1 = Molecule::new("mol1");
        mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

        let mol2 = mol1.clone();

        let rmsd = mol1.rmsd_to(&mol2).unwrap();
        assert!(rmsd.abs() < 1e-10);
    }

    #[test]
    fn test_molecule_apply_rotation_matrix_method() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        mol.apply_rotation_matrix(&identity);

        assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_molecule_apply_transform_method() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        mol.apply_transform(&identity, [1.0, 2.0, 3.0]);

        assert!((mol.atoms[0].x - 1.0).abs() < 1e-10);
        assert!((mol.atoms[0].y - 2.0).abs() < 1e-10);
        assert!((mol.atoms[0].z - 3.0).abs() < 1e-10);
    }
}
