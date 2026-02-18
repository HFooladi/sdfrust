//! Bond angle and dihedral angle calculations for 3D molecular structures.
//!
//! Provides angle computations needed by DimeNet, GemNet, and other
//! angular message-passing GNNs.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::geometry::angles;
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 0.9572, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -0.2400, 0.9266, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! let angle = angles::bond_angle(&mol, 1, 0, 2).unwrap();
//! // H-O-H angle should be ~104.5°
//! assert!(angle > 1.5 && angle < 2.0); // radians
//! ```

use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Result of enumerating all bond angles.
#[derive(Debug, Clone)]
pub struct AngleResult {
    /// Triplet indices: each row is (i, j, k) where j is the central atom.
    pub triplets: Vec<[usize; 3]>,
    /// Angle values in radians.
    pub angles: Vec<f64>,
}

/// Result of enumerating all dihedral angles.
#[derive(Debug, Clone)]
pub struct DihedralResult {
    /// Quadruplet indices: each row is (i, j, k, l).
    pub quadruplets: Vec<[usize; 4]>,
    /// Dihedral angle values in radians (-π to π).
    pub angles: Vec<f64>,
}

/// Compute the bond angle at atom j formed by atoms i-j-k (in radians).
///
/// Returns `None` if any atom index is out of bounds or if the angle
/// is undefined (e.g., degenerate geometry with zero-length vectors).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::angles;
///
/// let mut mol = Molecule::new("right_angle");
/// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
///
/// let angle = angles::bond_angle(&mol, 0, 1, 2).unwrap();
/// assert!((angle - std::f64::consts::FRAC_PI_2).abs() < 1e-10);
/// ```
pub fn bond_angle(mol: &Molecule, i: usize, j: usize, k: usize) -> Option<f64> {
    let ai = mol.atoms.get(i)?;
    let aj = mol.atoms.get(j)?;
    let ak = mol.atoms.get(k)?;

    // Vectors from j to i and j to k
    let v1 = [ai.x - aj.x, ai.y - aj.y, ai.z - aj.z];
    let v2 = [ak.x - aj.x, ak.y - aj.y, ak.z - aj.z];

    let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();

    let denom = len1 * len2;
    if denom < 1e-15 {
        return None;
    }

    let cos_angle = (dot / denom).clamp(-1.0, 1.0);
    Some(cos_angle.acos())
}

/// Compute the dihedral (torsion) angle for atoms i-j-k-l (in radians).
///
/// Returns a value in the range [-π, π]. Returns `None` if any atom
/// index is out of bounds or if the dihedral is undefined.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::angles;
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
/// mol.atoms.push(Atom::new(3, "C", 0.0, 1.0, 1.0));
///
/// let dihedral = angles::dihedral_angle(&mol, 0, 1, 2, 3);
/// assert!(dihedral.is_some());
/// ```
pub fn dihedral_angle(mol: &Molecule, i: usize, j: usize, k: usize, l: usize) -> Option<f64> {
    let ai = mol.atoms.get(i)?;
    let aj = mol.atoms.get(j)?;
    let ak = mol.atoms.get(k)?;
    let al = mol.atoms.get(l)?;

    // Vectors along the chain
    let b1 = [aj.x - ai.x, aj.y - ai.y, aj.z - ai.z];
    let b2 = [ak.x - aj.x, ak.y - aj.y, ak.z - aj.z];
    let b3 = [al.x - ak.x, al.y - ak.y, al.z - ak.z];

    // Normal vectors to the two planes
    let n1 = cross(&b1, &b2);
    let n2 = cross(&b2, &b3);

    let n1_len = vec_len(&n1);
    let n2_len = vec_len(&n2);

    if n1_len < 1e-15 || n2_len < 1e-15 {
        return None;
    }

    // Compute dihedral using atan2 for proper sign
    let m1 = cross(
        &n1,
        &[
            b2[0] / vec_len(&b2),
            b2[1] / vec_len(&b2),
            b2[2] / vec_len(&b2),
        ],
    );
    let x = dot(&n1, &n2);
    let y = dot(&m1, &n2);

    Some((-y).atan2(x))
}

/// Enumerate all bond angles from bonded paths.
///
/// For each atom j with at least 2 bonded neighbors, enumerates all
/// unique triplets (i, j, k) where both i-j and j-k are bonds.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::geometry::angles;
///
/// let mut mol = Molecule::new("water");
/// mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 0.9572, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -0.2400, 0.9266, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
///
/// let result = angles::all_bond_angles(&mol);
/// assert_eq!(result.triplets.len(), 1); // One H-O-H angle
/// ```
pub fn all_bond_angles(mol: &Molecule) -> AngleResult {
    let adj = AdjacencyList::from_molecule(mol);
    let mut triplets = Vec::new();
    let mut angles = Vec::new();

    for j in 0..adj.num_atoms() {
        let neighbors = adj.neighbor_atoms(j);
        for ni in 0..neighbors.len() {
            for nk in (ni + 1)..neighbors.len() {
                let i = neighbors[ni];
                let k = neighbors[nk];
                if let Some(angle) = bond_angle(mol, i, j, k) {
                    triplets.push([i, j, k]);
                    angles.push(angle);
                }
            }
        }
    }

    AngleResult { triplets, angles }
}

/// Enumerate all dihedral angles from bonded paths.
///
/// For each bond j-k, finds all quadruplets (i, j, k, l) where
/// i is bonded to j and l is bonded to k (excluding j-k bond itself).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::geometry::angles;
///
/// let mut mol = Molecule::new("butane_skel");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "C", 2.5, 1.2, 0.0));
/// mol.atoms.push(Atom::new(3, "C", 3.5, 1.2, 1.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
/// mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
///
/// let result = angles::all_dihedral_angles(&mol);
/// assert_eq!(result.quadruplets.len(), 1); // 0-1-2-3
/// ```
pub fn all_dihedral_angles(mol: &Molecule) -> DihedralResult {
    let adj = AdjacencyList::from_molecule(mol);
    let mut quadruplets = Vec::new();
    let mut angles = Vec::new();

    for bond in &mol.bonds {
        let j = bond.atom1;
        let k = bond.atom2;

        let neighbors_j = adj.neighbor_atoms(j);
        let neighbors_k = adj.neighbor_atoms(k);

        for &i in &neighbors_j {
            if i == k {
                continue;
            }
            for &l in &neighbors_k {
                if l == j || l == i {
                    continue;
                }
                if let Some(angle) = dihedral_angle(mol, i, j, k, l) {
                    quadruplets.push([i, j, k, l]);
                    angles.push(angle);
                }
            }
        }
    }

    DihedralResult {
        quadruplets,
        angles,
    }
}

/// Cross product of two 3D vectors.
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3D vectors.
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Length of a 3D vector.
fn vec_len(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn test_bond_angle_right_angle() {
        let mut mol = Molecule::new("right_angle");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));

        let angle = bond_angle(&mol, 0, 1, 2).unwrap();
        assert!((angle - FRAC_PI_2).abs() < 1e-10);
    }

    #[test]
    fn test_bond_angle_straight() {
        let mut mol = Molecule::new("linear");
        mol.atoms.push(Atom::new(0, "C", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 1.0, 0.0, 0.0));

        let angle = bond_angle(&mol, 0, 1, 2).unwrap();
        assert!((angle - PI).abs() < 1e-10);
    }

    #[test]
    fn test_bond_angle_out_of_bounds() {
        let mol = Molecule::new("empty");
        assert!(bond_angle(&mol, 0, 1, 2).is_none());
    }

    #[test]
    fn test_dihedral_angle_planar() {
        let mut mol = Molecule::new("planar");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 2.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 3.0, 1.0, 0.0));

        let dih = dihedral_angle(&mol, 0, 1, 2, 3).unwrap();
        // All atoms in the same plane → dihedral should be ~0 or ~π
        assert!(dih.abs() < 0.1 || (dih.abs() - PI).abs() < 0.1);
    }

    #[test]
    fn test_dihedral_angle_perpendicular() {
        let mut mol = Molecule::new("perp");
        mol.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 0.0, 1.0, 1.0));

        let dih = dihedral_angle(&mol, 0, 1, 2, 3).unwrap();
        assert!((dih.abs() - FRAC_PI_2).abs() < 0.1);
    }

    #[test]
    fn test_all_bond_angles_water() {
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.9572, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.2400, 0.9266, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

        let result = all_bond_angles(&mol);
        assert_eq!(result.triplets.len(), 1);
        // Should be the H-O-H angle
        assert_eq!(result.triplets[0][1], 0); // O is center
    }

    #[test]
    fn test_all_dihedral_angles_butane() {
        let mut mol = Molecule::new("butane_skel");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 2.5, 1.2, 0.0));
        mol.atoms.push(Atom::new(3, "C", 3.5, 1.2, 1.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));

        let result = all_dihedral_angles(&mol);
        assert_eq!(result.quadruplets.len(), 1);
    }

    #[test]
    fn test_all_bond_angles_methane() {
        let mut mol = Molecule::new("methane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.6289, 0.6289, 0.6289));
        mol.atoms.push(Atom::new(2, "H", -0.6289, -0.6289, 0.6289));
        mol.atoms.push(Atom::new(3, "H", -0.6289, 0.6289, -0.6289));
        mol.atoms.push(Atom::new(4, "H", 0.6289, -0.6289, -0.6289));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let result = all_bond_angles(&mol);
        // C(4) = 6 unique pairs of neighbors
        assert_eq!(result.triplets.len(), 6);
        // All angles should be ~109.47° (tetrahedral)
        let expected = 109.47_f64.to_radians();
        for angle in &result.angles {
            assert!(
                (angle - expected).abs() < 0.02,
                "Expected ~{:.2}°, got {:.2}°",
                expected.to_degrees(),
                angle.to_degrees()
            );
        }
    }

    #[test]
    fn test_empty_molecule_angles() {
        let mol = Molecule::new("empty");
        let result = all_bond_angles(&mol);
        assert!(result.triplets.is_empty());

        let result = all_dihedral_angles(&mol);
        assert!(result.quadruplets.is_empty());
    }
}
