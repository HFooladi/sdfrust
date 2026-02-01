//! RMSD (Root Mean Square Deviation) calculation.

use crate::{Molecule, SdfError};

/// Calculate RMSD between two sets of coordinates.
///
/// The coordinates are provided as slices of (x, y, z) tuples.
/// Both coordinate sets must have the same length.
///
/// # Arguments
///
/// * `coords1` - First set of coordinates
/// * `coords2` - Second set of coordinates
///
/// # Returns
///
/// The RMSD value, or an error if the coordinate sets have different lengths.
///
/// # Example
///
/// ```rust
/// use sdfrust::geometry::rmsd_from_coords;
///
/// let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
/// let coords2 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
/// let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
/// assert!(rmsd.abs() < 1e-10);
/// ```
pub fn rmsd_from_coords(
    coords1: &[(f64, f64, f64)],
    coords2: &[(f64, f64, f64)],
) -> Result<f64, SdfError> {
    if coords1.len() != coords2.len() {
        return Err(SdfError::AtomCountMismatch {
            expected: coords1.len(),
            found: coords2.len(),
        });
    }

    if coords1.is_empty() {
        return Ok(0.0);
    }

    let n = coords1.len() as f64;
    let sum_sq: f64 = coords1
        .iter()
        .zip(coords2.iter())
        .map(|((x1, y1, z1), (x2, y2, z2))| {
            let dx = x1 - x2;
            let dy = y1 - y2;
            let dz = z1 - z2;
            dx * dx + dy * dy + dz * dz
        })
        .sum();

    Ok((sum_sq / n).sqrt())
}

/// Calculate RMSD between two molecules.
///
/// Computes the root mean square deviation of atomic positions
/// between two molecules. The molecules must have the same number of atoms.
/// No alignment is performed - atoms are compared directly by index.
///
/// # Arguments
///
/// * `mol1` - First molecule
/// * `mol2` - Second molecule
///
/// # Returns
///
/// The RMSD value in Angstroms, or an error if the molecules have different atom counts.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::rmsd_to;
///
/// let mut mol1 = Molecule::new("mol1");
/// mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol1.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
///
/// let mut mol2 = Molecule::new("mol2");
/// mol2.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol2.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
///
/// let rmsd = rmsd_to(&mol1, &mol2).unwrap();
/// assert!(rmsd.abs() < 1e-10);
/// ```
pub fn rmsd_to(mol1: &Molecule, mol2: &Molecule) -> Result<f64, SdfError> {
    let coords1: Vec<(f64, f64, f64)> = mol1.atoms.iter().map(|a| (a.x, a.y, a.z)).collect();
    let coords2: Vec<(f64, f64, f64)> = mol2.atoms.iter().map(|a| (a.x, a.y, a.z)).collect();
    rmsd_from_coords(&coords1, &coords2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;

    #[test]
    fn test_rmsd_identical_molecules() {
        let mut mol1 = Molecule::new("mol1");
        mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol1.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol1.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));

        let mol2 = mol1.clone();

        let rmsd = rmsd_to(&mol1, &mol2).unwrap();
        assert!(rmsd.abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_translated_molecule() {
        let mut mol1 = Molecule::new("mol1");
        mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol1.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));

        let mut mol2 = Molecule::new("mol2");
        mol2.atoms.push(Atom::new(0, "C", 1.0, 0.0, 0.0)); // translated by (1, 0, 0)
        mol2.atoms.push(Atom::new(1, "C", 2.0, 0.0, 0.0));

        let rmsd = rmsd_to(&mol1, &mol2).unwrap();
        // Each atom is displaced by 1.0, so RMSD = sqrt(1^2 + 1^2) / sqrt(2) = 1.0
        assert!((rmsd - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_different_atom_counts() {
        let mut mol1 = Molecule::new("mol1");
        mol1.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));

        let mut mol2 = Molecule::new("mol2");
        mol2.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol2.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));

        let result = rmsd_to(&mol1, &mol2);
        assert!(result.is_err());
    }

    #[test]
    fn test_rmsd_empty_molecules() {
        let mol1 = Molecule::new("mol1");
        let mol2 = Molecule::new("mol2");

        let rmsd = rmsd_to(&mol1, &mol2).unwrap();
        assert!(rmsd.abs() < 1e-10);
    }

    #[test]
    fn test_rmsd_from_coords() {
        let coords1 = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)];
        let coords2 = vec![(0.5, 0.0, 0.0), (1.5, 0.0, 0.0)];

        let rmsd = rmsd_from_coords(&coords1, &coords2).unwrap();
        // Each atom displaced by 0.5
        assert!((rmsd - 0.5).abs() < 1e-10);
    }
}
