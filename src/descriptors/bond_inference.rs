//! Bond inference from 3D coordinates.
//!
//! Infers connectivity (single bonds) from atomic positions and covalent radii,
//! similar to the approach used by [xyz2mol](https://github.com/jensengroup/xyz2mol)
//! and Open Babel. Two atoms are considered bonded if their distance is within
//! the sum of their covalent radii plus a tolerance.
//!
//! This is particularly useful for XYZ files, which contain only coordinates
//! and element types without bond information.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, infer_bonds};
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
//! mol.atoms.push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
//! mol.atoms.push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));
//!
//! infer_bonds(&mut mol, None).unwrap();
//! assert_eq!(mol.bond_count(), 2);
//! ```

use crate::bond::{Bond, BondOrder};
use crate::descriptors::elements::covalent_radius;
use crate::error::{Result, SdfError};
use crate::molecule::Molecule;

/// Default tolerance in Angstroms added to the sum of covalent radii.
///
/// This value (0.45 A) matches the default used by xyz2mol and Open Babel.
pub const DEFAULT_TOLERANCE: f64 = 0.45;

/// Minimum distance threshold to avoid bonding overlapping atoms.
const MIN_DISTANCE: f64 = 0.01;

/// Configuration for bond inference.
#[derive(Debug, Clone)]
pub struct BondInferenceConfig {
    /// Tolerance in Angstroms added to the sum of covalent radii (default: 0.45).
    pub tolerance: f64,
    /// Whether to clear existing bonds before inference (default: true).
    pub clear_existing_bonds: bool,
}

impl Default for BondInferenceConfig {
    fn default() -> Self {
        Self {
            tolerance: DEFAULT_TOLERANCE,
            clear_existing_bonds: true,
        }
    }
}

/// Infer single bonds from 3D coordinates and covalent radii.
///
/// Two atoms are bonded if their distance is within the sum of their
/// covalent radii plus a tolerance. All inferred bonds are single bonds.
///
/// # Arguments
///
/// * `mol` - The molecule to add bonds to (existing bonds are cleared by default)
/// * `tolerance` - Optional tolerance in Angstroms (default: 0.45)
///
/// # Errors
///
/// Returns `SdfError::BondInferenceError` if any atom has an unknown element
/// with no covalent radius data.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, infer_bonds};
///
/// let mut mol = Molecule::new("water");
/// mol.atoms.push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
/// mol.atoms.push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
/// mol.atoms.push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));
///
/// infer_bonds(&mut mol, None).unwrap();
/// assert_eq!(mol.bond_count(), 2);
/// ```
pub fn infer_bonds(mol: &mut Molecule, tolerance: Option<f64>) -> Result<()> {
    let config = BondInferenceConfig {
        tolerance: tolerance.unwrap_or(DEFAULT_TOLERANCE),
        ..Default::default()
    };
    infer_bonds_with_config(mol, &config)
}

/// Infer single bonds from 3D coordinates with full configuration.
///
/// See [`infer_bonds`] for a simpler API.
///
/// # Arguments
///
/// * `mol` - The molecule to add bonds to
/// * `config` - Configuration controlling tolerance and behavior
///
/// # Errors
///
/// Returns `SdfError::BondInferenceError` if any atom has an unknown element
/// with no covalent radius data.
pub fn infer_bonds_with_config(mol: &mut Molecule, config: &BondInferenceConfig) -> Result<()> {
    // Look up covalent radii for all atoms
    let radii: Vec<f64> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            covalent_radius(&atom.element).ok_or_else(|| SdfError::BondInferenceError {
                element: atom.element.clone(),
                index: i,
            })
        })
        .collect::<Result<Vec<f64>>>()?;

    if config.clear_existing_bonds {
        mol.bonds.clear();
    }

    let n = mol.atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = mol.atoms[i].x - mol.atoms[j].x;
            let dy = mol.atoms[i].y - mol.atoms[j].y;
            let dz = mol.atoms[i].z - mol.atoms[j].z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            let max_dist = radii[i] + radii[j] + config.tolerance;

            if dist > MIN_DISTANCE && dist <= max_dist {
                mol.bonds.push(Bond::new(i, j, BondOrder::Single));
            }
        }
    }

    Ok(())
}

// Extension methods for Molecule
impl Molecule {
    /// Infer single bonds from 3D coordinates and covalent radii.
    ///
    /// Two atoms are bonded if their distance is within the sum of their
    /// covalent radii plus a tolerance. All inferred bonds are single bonds.
    /// Existing bonds are cleared before inference.
    ///
    /// # Arguments
    ///
    /// * `tolerance` - Optional tolerance in Angstroms (default: 0.45)
    ///
    /// # Errors
    ///
    /// Returns `SdfError::BondInferenceError` if any atom has an unknown element.
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom};
    ///
    /// let mut mol = Molecule::new("water");
    /// mol.atoms.push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
    /// mol.atoms.push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
    /// mol.atoms.push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));
    ///
    /// mol.infer_bonds(None).unwrap();
    /// assert_eq!(mol.bond_count(), 2);
    /// ```
    pub fn infer_bonds(&mut self, tolerance: Option<f64>) -> Result<()> {
        infer_bonds(self, tolerance)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;

    #[test]
    fn test_infer_bonds_water() {
        let mut mol = Molecule::new("water");
        mol.atoms
            .push(Atom::new(0, "O", 0.000000, 0.000000, 0.117300));
        mol.atoms
            .push(Atom::new(1, "H", 0.756950, 0.000000, -0.469200));
        mol.atoms
            .push(Atom::new(2, "H", -0.756950, 0.000000, -0.469200));

        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 2);
        assert!(mol.bonds.iter().all(|b| b.order == BondOrder::Single));
    }

    #[test]
    fn test_infer_bonds_empty_molecule() {
        let mut mol = Molecule::new("empty");
        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn test_infer_bonds_single_atom() {
        let mut mol = Molecule::new("single");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn test_infer_bonds_unknown_element() {
        let mut mol = Molecule::new("unknown");
        mol.atoms.push(Atom::new(0, "Xx", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        let result = infer_bonds(&mut mol, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_infer_bonds_distant_atoms() {
        let mut mol = Molecule::new("distant");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 10.0, 0.0, 0.0));
        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 0);
    }

    #[test]
    fn test_infer_bonds_clears_existing() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Single);
    }

    #[test]
    fn test_infer_bonds_keep_existing() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        let config = BondInferenceConfig {
            clear_existing_bonds: false,
            ..Default::default()
        };
        infer_bonds_with_config(&mut mol, &config).unwrap();
        // Original double bond + inferred single bond
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn test_infer_bonds_tolerance_effect() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.9, 0.0, 0.0));

        // With default tolerance (0.45), C-C sum = 0.77+0.77+0.45 = 1.99 → bonded
        infer_bonds(&mut mol, None).unwrap();
        assert_eq!(mol.bond_count(), 1);

        // With very small tolerance, not bonded
        infer_bonds(&mut mol, Some(0.1)).unwrap();
        assert_eq!(mol.bond_count(), 0);
    }
}
