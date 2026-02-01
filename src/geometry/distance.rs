//! Distance calculations for molecular structures.

use crate::Molecule;

/// Compute the pairwise distance matrix for all atoms in a molecule.
///
/// Returns an NxN matrix where entry \[i\]\[j\] is the Euclidean distance
/// between atom i and atom j in Angstroms.
///
/// # Arguments
///
/// * `molecule` - The molecule to compute distances for
///
/// # Returns
///
/// A symmetric matrix of pairwise distances. The diagonal is all zeros.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::distance_matrix;
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
///
/// let matrix = distance_matrix(&mol);
/// assert_eq!(matrix.len(), 3);
/// assert!((matrix[0][1] - 1.0).abs() < 1e-10);
/// ```
pub fn distance_matrix(molecule: &Molecule) -> Vec<Vec<f64>> {
    let n = molecule.atoms.len();
    let mut matrix = vec![vec![0.0; n]; n];

    for (i, atom_i) in molecule.atoms.iter().enumerate() {
        for (j, atom_j) in molecule.atoms.iter().enumerate().skip(i + 1) {
            let dist = atom_i.distance_to(atom_j);
            matrix[i][j] = dist;
            matrix[j][i] = dist;
        }
    }

    matrix
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;

    #[test]
    fn test_distance_matrix_empty() {
        let mol = Molecule::new("empty");
        let matrix = distance_matrix(&mol);
        assert!(matrix.is_empty());
    }

    #[test]
    fn test_distance_matrix_single_atom() {
        let mut mol = Molecule::new("single");
        mol.atoms.push(Atom::new(0, "C", 1.0, 2.0, 3.0));

        let matrix = distance_matrix(&mol);
        assert_eq!(matrix.len(), 1);
        assert_eq!(matrix[0].len(), 1);
        assert!((matrix[0][0]).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_two_atoms() {
        let mut mol = Molecule::new("two");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 3.0, 4.0, 0.0));

        let matrix = distance_matrix(&mol);
        assert_eq!(matrix.len(), 2);
        assert!((matrix[0][0]).abs() < 1e-10);
        assert!((matrix[1][1]).abs() < 1e-10);
        assert!((matrix[0][1] - 5.0).abs() < 1e-10);
        assert!((matrix[1][0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance_matrix_symmetric() {
        let mut mol = Molecule::new("three");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));

        let matrix = distance_matrix(&mol);
        assert_eq!(matrix.len(), 3);

        // Check symmetry
        for (i, row) in matrix.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - matrix[j][i]).abs() < 1e-10);
            }
        }
    }
}
