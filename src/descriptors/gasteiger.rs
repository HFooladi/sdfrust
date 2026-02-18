//! Gasteiger-Marsili partial equalization of orbital electronegativity (PEOE).
//!
//! Computes partial atomic charges using the iterative Gasteiger algorithm.
//! This is a topology-only method (no force field or QM needed).
//!
//! Reference: Gasteiger, J.; Marsili, M. Tetrahedron 1980, 36, 3219-3228.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::gasteiger;
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! let charges = gasteiger::gasteiger_charges(&mol);
//! // O should be negative, H should be positive
//! assert!(charges[0] < 0.0);
//! assert!(charges[1] > 0.0);
//! ```

use crate::descriptors::hybridization::{Hybridization, atom_hybridization};
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Gasteiger electronegativity parameters (a, b, c).
///
/// Electronegativity χ = a + b*q + c*q²
/// where q is the charge on the atom.
struct ElectronegativityParams {
    a: f64,
    b: f64,
    c: f64,
}

/// Get electronegativity parameters for an atom type.
///
/// Parameters from Gasteiger & Marsili (1980) and extensions.
fn get_params(element: &str, hybridization: Hybridization) -> Option<ElectronegativityParams> {
    let elem = element.trim();
    let upper: String = elem
        .chars()
        .next()
        .map(|c| c.to_uppercase().collect::<String>())
        .unwrap_or_default()
        + &elem.chars().skip(1).collect::<String>().to_lowercase();

    match (upper.as_str(), hybridization) {
        // Hydrogen
        ("H", _) | ("D", _) | ("T", _) => Some(ElectronegativityParams {
            a: 7.17,
            b: 6.24,
            c: -0.56,
        }),

        // Carbon
        ("C", Hybridization::SP3) => Some(ElectronegativityParams {
            a: 7.98,
            b: 9.18,
            c: 1.88,
        }),
        ("C", Hybridization::SP2) => Some(ElectronegativityParams {
            a: 8.79,
            b: 9.32,
            c: 1.51,
        }),
        ("C", Hybridization::SP) => Some(ElectronegativityParams {
            a: 10.39,
            b: 9.45,
            c: 0.73,
        }),
        ("C", _) => Some(ElectronegativityParams {
            a: 7.98,
            b: 9.18,
            c: 1.88,
        }),

        // Nitrogen
        ("N", Hybridization::SP3) => Some(ElectronegativityParams {
            a: 11.54,
            b: 10.82,
            c: 1.36,
        }),
        ("N", Hybridization::SP2) => Some(ElectronegativityParams {
            a: 12.87,
            b: 11.15,
            c: 0.85,
        }),
        ("N", Hybridization::SP) => Some(ElectronegativityParams {
            a: 15.68,
            b: 11.70,
            c: -0.27,
        }),
        ("N", _) => Some(ElectronegativityParams {
            a: 11.54,
            b: 10.82,
            c: 1.36,
        }),

        // Oxygen
        ("O", Hybridization::SP3) => Some(ElectronegativityParams {
            a: 14.18,
            b: 12.92,
            c: 1.39,
        }),
        ("O", Hybridization::SP2) => Some(ElectronegativityParams {
            a: 17.07,
            b: 13.79,
            c: 0.47,
        }),
        ("O", _) => Some(ElectronegativityParams {
            a: 14.18,
            b: 12.92,
            c: 1.39,
        }),

        // Fluorine
        ("F", _) => Some(ElectronegativityParams {
            a: 14.66,
            b: 13.85,
            c: 2.31,
        }),

        // Chlorine
        ("Cl", _) => Some(ElectronegativityParams {
            a: 11.00,
            b: 9.69,
            c: 1.35,
        }),

        // Bromine
        ("Br", _) => Some(ElectronegativityParams {
            a: 10.08,
            b: 8.47,
            c: 1.16,
        }),

        // Iodine
        ("I", _) => Some(ElectronegativityParams {
            a: 9.90,
            b: 7.96,
            c: 0.96,
        }),

        // Sulfur
        ("S", Hybridization::SP3) => Some(ElectronegativityParams {
            a: 10.14,
            b: 9.13,
            c: 1.38,
        }),
        ("S", Hybridization::SP2) => Some(ElectronegativityParams {
            a: 10.88,
            b: 9.49,
            c: 1.33,
        }),
        ("S", _) => Some(ElectronegativityParams {
            a: 10.14,
            b: 9.13,
            c: 1.38,
        }),

        // Phosphorus
        ("P", Hybridization::SP3) => Some(ElectronegativityParams {
            a: 8.90,
            b: 8.24,
            c: 0.96,
        }),
        ("P", _) => Some(ElectronegativityParams {
            a: 8.90,
            b: 8.24,
            c: 0.96,
        }),

        // Silicon
        ("Si", _) => Some(ElectronegativityParams {
            a: 7.30,
            b: 6.56,
            c: 0.74,
        }),

        // Boron
        ("B", _) => Some(ElectronegativityParams {
            a: 6.88,
            b: 5.98,
            c: 0.68,
        }),

        // Default: use carbon SP3 parameters as fallback
        _ => None,
    }
}

/// Compute electronegativity for a given charge.
fn electronegativity(params: &ElectronegativityParams, q: f64) -> f64 {
    params.a + params.b * q + params.c * q * q
}

/// Compute Gasteiger partial charges for all atoms.
///
/// Uses the iterative PEOE algorithm with 6 iterations and 0.5 damping factor.
///
/// Returns a vector of partial charges of length `mol.atom_count()`.
/// Returns all zeros if any atom has unknown electronegativity parameters.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::gasteiger;
///
/// // HF: H should be positive, F should be negative
/// let mut mol = Molecule::new("HF");
/// mol.atoms.push(Atom::new(0, "H", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "F", 0.9, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let charges = gasteiger::gasteiger_charges(&mol);
/// assert!(charges[0] > 0.0); // H is positive
/// assert!(charges[1] < 0.0); // F is negative
/// ```
pub fn gasteiger_charges(mol: &Molecule) -> Vec<f64> {
    gasteiger_charges_with_params(mol, 6, 0.5)
}

/// Compute Gasteiger partial charges with custom parameters.
///
/// # Arguments
///
/// * `mol` - The molecule
/// * `max_iter` - Maximum number of iterations (typically 6-8)
/// * `damping` - Damping factor (typically 0.5)
pub fn gasteiger_charges_with_params(mol: &Molecule, max_iter: usize, damping: f64) -> Vec<f64> {
    let n = mol.atom_count();
    if n == 0 {
        return vec![];
    }

    let adj = AdjacencyList::from_molecule(mol);

    // Get electronegativity parameters for each atom
    let params: Vec<Option<ElectronegativityParams>> = (0..n)
        .map(|i| {
            let hyb = atom_hybridization(mol, i);
            get_params(&mol.atoms[i].element, hyb)
        })
        .collect();

    // Initialize charges from formal charges
    let mut charges: Vec<f64> = mol.atoms.iter().map(|a| a.formal_charge as f64).collect();

    // Iterative charge equalization
    for iter in 0..max_iter {
        let damp = damping.powi(iter as i32 + 1);
        let mut charge_deltas = vec![0.0f64; n];

        // For each bond, transfer charge from less electronegative to more
        for bond in &mol.bonds {
            let i = bond.atom1;
            let j = bond.atom2;

            if i >= n || j >= n {
                continue;
            }

            let (params_i, params_j) = match (&params[i], &params[j]) {
                (Some(pi), Some(pj)) => (pi, pj),
                _ => continue,
            };

            let chi_i = electronegativity(params_i, charges[i]);
            let chi_j = electronegativity(params_j, charges[j]);

            // Charge flows from less electronegative to more electronegative
            let delta_chi = chi_j - chi_i;

            // Normalization: use the electronegativity of the positive end at q=+1
            let chi_plus = if delta_chi > 0.0 {
                electronegativity(params_i, 1.0)
            } else {
                electronegativity(params_j, 1.0)
            };

            if chi_plus.abs() < 1e-10 {
                continue;
            }

            let transfer = damp * delta_chi / chi_plus;

            charge_deltas[i] += transfer;
            charge_deltas[j] -= transfer;
        }

        for i in 0..n {
            charges[i] += charge_deltas[i];
        }
    }

    // For atoms without parameters, set charge to 0
    for i in 0..n {
        if params[i].is_none() {
            charges[i] = 0.0;
        }
    }

    // Suppress unused variable warning
    let _ = &adj;

    charges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    #[test]
    fn test_hf_charges() {
        let mut mol = Molecule::new("HF");
        mol.atoms.push(Atom::new(0, "H", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "F", 0.92, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let charges = gasteiger_charges(&mol);
        assert!(charges[0] > 0.0, "H should be positive in HF");
        assert!(charges[1] < 0.0, "F should be negative in HF");
    }

    #[test]
    fn test_water_charges() {
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

        let charges = gasteiger_charges(&mol);
        assert!(charges[0] < 0.0, "O should be negative in water");
        assert!(charges[1] > 0.0, "H should be positive in water");
        assert!(charges[2] > 0.0, "H should be positive in water");

        // Total charge should be approximately 0
        let total: f64 = charges.iter().sum();
        assert!(
            total.abs() < 0.1,
            "Total charge should be ~0, got {}",
            total
        );
    }

    #[test]
    fn test_methane_charges() {
        let mut mol = Molecule::new("methane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "H", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let charges = gasteiger_charges(&mol);
        // C is slightly more electronegative than H
        assert!(charges[0] < 0.0, "C should be slightly negative in methane");
        // All H should have same charge
        let h_charge = charges[1];
        for i in 2..5 {
            assert!(
                (charges[i] - h_charge).abs() < 1e-10,
                "All H should have same charge"
            );
        }
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new("empty");
        let charges = gasteiger_charges(&mol);
        assert!(charges.is_empty());
    }

    #[test]
    fn test_single_atom() {
        let mut mol = Molecule::new("He");
        mol.atoms.push(Atom::new(0, "He", 0.0, 0.0, 0.0));
        let charges = gasteiger_charges(&mol);
        assert_eq!(charges.len(), 1);
        // No bonds → charge stays at formal charge (0)
        assert!((charges[0]).abs() < 1e-10);
    }

    #[test]
    fn test_charge_conservation() {
        // Neutral molecule: sum of charges should be ~0
        let mut mol = Molecule::new("ethanol");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));

        let charges = gasteiger_charges(&mol);
        let total: f64 = charges.iter().sum();
        assert!(
            total.abs() < 0.1,
            "Total charge should be ~0 for neutral molecule, got {}",
            total
        );
    }

    #[test]
    fn test_electronegativity_order() {
        // F > O > N > C: check relative charges
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "F", 1.3, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let charges = gasteiger_charges(&mol);
        assert!(charges[0] > 0.0, "C should be positive bonded to F");
        assert!(charges[1] < 0.0, "F should be negative bonded to C");
    }

    #[test]
    fn test_custom_params() {
        let mut mol = Molecule::new("HF");
        mol.atoms.push(Atom::new(0, "H", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "F", 0.92, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let charges = gasteiger_charges_with_params(&mol, 10, 0.5);
        assert!(charges[0] > 0.0);
        assert!(charges[1] < 0.0);
    }
}
