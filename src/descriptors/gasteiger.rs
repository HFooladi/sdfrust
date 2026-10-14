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

use crate::descriptors::conjugation::all_conjugated_bonds;
use crate::descriptors::elements::get_element;
use crate::descriptors::hybridization::{Hybridization, all_hybridizations};
use crate::descriptors::valence::all_implicit_hydrogen_counts;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Damping factor applied (and squared each iteration) to charge transfers.
const DAMPING: f64 = 0.5;

/// Ionization electronegativity of hydrogen (χ at q = +1), a special value
/// in the Gasteiger-Marsili scheme.
const HYDROGEN_ION_X: f64 = 20.02;

/// Default number of iterations, as in RDKit.
const DEFAULT_ITERATIONS: usize = 12;

/// Gasteiger electronegativity parameters: χ(q) = a + b·q + c·q².
#[derive(Clone, Copy)]
struct ElectronegativityParams {
    a: f64,
    b: f64,
    c: f64,
}

/// Parameter table (element, mode, a, b, c), identical to RDKit's
/// `GasteigerParams` defaults. Modes are the hybridization ("sp3", "sp2",
/// "sp"), "*" for hydrogen, and "so"/"so2" for sulfur without a hybridization.
const PARAMS: &[(&str, &str, f64, f64, f64)] = &[
    ("H", "*", 7.17, 6.24, -0.56),
    ("C", "sp3", 7.98, 9.18, 1.88),
    ("C", "sp2", 8.79, 9.32, 1.51),
    ("C", "sp", 10.39, 9.45, 0.73),
    ("N", "sp3", 11.54, 10.82, 1.36),
    ("N", "sp2", 12.87, 11.15, 0.85),
    ("N", "sp", 15.68, 11.7, -0.27),
    ("O", "sp3", 14.18, 12.92, 1.39),
    ("O", "sp2", 17.07, 13.79, 0.47),
    ("F", "sp3", 14.66, 13.85, 2.31),
    ("Cl", "sp3", 11.00, 9.69, 1.35),
    ("Br", "sp3", 10.08, 8.47, 1.16),
    ("I", "sp3", 9.9, 7.96, 0.96),
    ("S", "sp3", 10.14, 9.13, 1.38),
    ("S", "so", 10.14, 9.13, 1.38),
    ("S", "so2", 12.00, 10.81, 1.20),
    ("S", "sp2", 10.88, 9.49, 1.33),
    ("P", "sp3", 8.90, 8.24, 0.96),
    ("P", "sp2", 9.665, 8.530, 0.735),
    ("Si", "sp3", 7.300, 6.567, 0.657),
    ("Si", "sp2", 7.905, 6.748, 0.443),
    ("Si", "sp", 9.065, 7.027, -0.002),
    ("B", "sp3", 5.980, 6.820, 1.605),
    ("B", "sp2", 6.420, 6.807, 1.322),
    ("Be", "sp3", 3.845, 6.755, 3.165),
    ("Be", "sp2", 4.005, 6.725, 3.035),
    ("Mg", "sp2", 3.565, 5.572, 2.197),
    ("Mg", "sp3", 3.300, 5.587, 2.447),
    ("Mg", "sp", 4.040, 5.472, 1.823),
    ("Al", "sp3", 5.375, 4.953, 0.867),
    ("Al", "sp2", 5.795, 5.020, 0.695),
];

/// Looks up the parameters for an element/mode pair.
fn lookup(symbol: &str, mode: &str) -> Option<ElectronegativityParams> {
    PARAMS
        .iter()
        .find(|(e, m, ..)| *e == symbol && *m == mode)
        .map(|&(_, _, a, b, c)| ElectronegativityParams { a, b, c })
}

/// Parameters for atom `idx`, choosing the mode the way RDKit does.
fn atom_params(
    mol: &Molecule,
    adj: &AdjacencyList,
    idx: usize,
    hybridization: Hybridization,
) -> Option<ElectronegativityParams> {
    let element = get_element(&mol.atoms[idx].element)?;
    let mode = match hybridization {
        Hybridization::SP3 => "sp3",
        Hybridization::SP2 => "sp2",
        Hybridization::SP => "sp",
        _ if element.atomic_number == 1 => "*",
        _ if element.atomic_number == 16 => {
            let oxygens = adj
                .neighbor_atoms(idx)
                .iter()
                .filter(|&&n| {
                    get_element(&mol.atoms[n].element).is_some_and(|e| e.atomic_number == 8)
                })
                .count();
            match oxygens {
                2 => "so2",
                1 => "so",
                _ => "sp3",
            }
        }
        _ => return None,
    };
    // Hydrogen isotopes share hydrogen's parameters
    let symbol = if element.atomic_number == 1 {
        "H"
    } else {
        element.symbol
    };
    lookup(
        symbol,
        if element.atomic_number == 1 {
            "*"
        } else {
            mode
        },
    )
}

/// Initial charges: formal charges, with the charge of an atom spread evenly
/// over same-element atoms two conjugated bonds away (e.g. both carboxylate
/// oxygens get -0.5), as in RDKit's `splitChargeConjugated`.
fn initial_charges(mol: &Molecule, adj: &AdjacencyList, conjugated: &[bool]) -> Vec<f64> {
    let n = mol.atom_count();
    let mut charges = vec![0.0; n];
    for i in 0..n {
        let atom = &mol.atoms[i];
        if atom.formal_charge == 0 || charges[i] != 0.0 {
            continue;
        }
        let mut total = atom.formal_charge as f64;
        let mut members = vec![i];
        for &(middle, b1) in adj.neighbors(i) {
            if !conjugated[b1] {
                continue;
            }
            for &(other, b2) in adj.neighbors(middle) {
                if b2 != b1 && conjugated[b2] && mol.atoms[other].element == atom.element {
                    total += mol.atoms[other].formal_charge as f64;
                    members.push(other);
                }
            }
        }
        let share = total / members.len() as f64;
        for m in members {
            charges[m] = share;
        }
    }
    charges
}

/// Compute Gasteiger partial charges for all atoms.
///
/// Port of RDKit's `ComputeGasteigerCharges` (12 iterations): implicit
/// hydrogens take part in the charge equalization, hydrogen uses the special
/// ionization electronegativity 20.02, parameters depend on the perceived
/// hybridization, and formal charges are first spread over conjugated
/// same-element atoms. The returned charge of each atom excludes its implicit
/// hydrogens (like RDKit's `_GasteigerCharge`); see
/// [`gasteiger_charges_with_hydrogens`] for those.
///
/// Atoms without parameters (most metals) keep their formal charge and do not
/// exchange charge; RDKit reports NaN for them instead.
///
/// Returns a vector of partial charges of length `mol.atom_count()`.
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
    gasteiger_charges_with_params(mol, DEFAULT_ITERATIONS, DAMPING)
}

/// Gasteiger charges plus the total charge carried by each atom's implicit
/// hydrogens: returns `(atom_charges, implicit_h_charges)`.
///
/// Summing both vectors gives the molecule's net formal charge.
pub fn gasteiger_charges_with_hydrogens(mol: &Molecule) -> (Vec<f64>, Vec<f64>) {
    compute(mol, DEFAULT_ITERATIONS, DAMPING)
}

/// Compute Gasteiger partial charges with custom parameters.
///
/// # Arguments
///
/// * `mol` - The molecule
/// * `max_iter` - Number of iterations (RDKit's default is 12)
/// * `damping` - Damping factor, squared each iteration (typically 0.5)
pub fn gasteiger_charges_with_params(mol: &Molecule, max_iter: usize, damping: f64) -> Vec<f64> {
    compute(mol, max_iter, damping).0
}

fn compute(mol: &Molecule, max_iter: usize, damping: f64) -> (Vec<f64>, Vec<f64>) {
    let n = mol.atom_count();
    if n == 0 {
        return (vec![], vec![]);
    }

    let adj = AdjacencyList::from_molecule(mol);
    let hybridizations = all_hybridizations(mol);
    let conjugated = all_conjugated_bonds(mol);
    let implicit_h = all_implicit_hydrogen_counts(mol);
    let params: Vec<Option<ElectronegativityParams>> = (0..n)
        .map(|i| atom_params(mol, &adj, i, hybridizations[i]))
        .collect();
    let hydrogen = lookup("H", "*").expect("hydrogen parameters");
    let ion_x: Vec<f64> = (0..n)
        .map(|i| match &params[i] {
            _ if get_element(&mol.atoms[i].element).is_some_and(|e| e.atomic_number == 1) => {
                HYDROGEN_ION_X
            }
            Some(p) => p.a + p.b + p.c,
            None => 0.0,
        })
        .collect();

    let mut charges = initial_charges(mol, &adj, &conjugated);
    let mut h_charges = vec![0.0; n];
    let mut energy = vec![0.0; n];
    let mut damp = damping;

    for _ in 0..max_iter {
        for i in 0..n {
            if let Some(p) = &params[i] {
                energy[i] = p.a + charges[i] * (p.b + p.c * charges[i]);
            }
        }
        for i in 0..n {
            if params[i].is_none() {
                continue;
            }
            let mut dq = 0.0;
            for &(j, _) in adj.neighbors(i) {
                if params[j].is_none() {
                    continue;
                }
                // Normalize by χ+ of the electron donor (the less electronegative atom)
                let dx = energy[j] - energy[i];
                dq += dx / if dx < 0.0 { ion_x[j] } else { ion_x[i] };
            }
            // Implicit hydrogens exchange charge with their heavy atom only,
            // so their charge is updated at the same time.
            let num_h = implicit_h[i] as f64;
            if num_h > 0.0 {
                let q_h = h_charges[i] / num_h;
                let h_energy = hydrogen.a + q_h * (hydrogen.b + hydrogen.c * q_h);
                let dx = h_energy - energy[i];
                let dq_h = dx / if dx < 0.0 { HYDROGEN_ION_X } else { ion_x[i] };
                dq += num_h * dq_h;
                h_charges[i] -= num_h * dq_h * damp;
            }
            charges[i] += damp * dq;
        }
        damp *= damping;
    }

    (charges, h_charges)
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
        for charge in charges.iter().take(5).skip(2) {
            assert!(
                (charge - h_charge).abs() < 1e-10,
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

        // Heavy-atom-only ethanol: implicit hydrogens carry part of the charge
        let (charges, h_charges) = gasteiger_charges_with_hydrogens(&mol);
        let total: f64 = charges.iter().chain(h_charges.iter()).sum();
        assert!(
            total.abs() < 1e-9,
            "Total charge should be 0 for neutral molecule, got {}",
            total
        );
        assert!(h_charges.iter().all(|&q| q > 0.0));
    }

    #[test]
    fn test_conjugated_formal_charge_is_split() {
        // Acetate: both oxygens start at -0.5 and end up equal
        let mut mol = Molecule::new("acetate");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 2.2, 1.0, 0.0));
        let mut o = Atom::new(3, "O", 2.2, -1.0, 0.0);
        o.formal_charge = -1;
        mol.atoms.push(o);
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
        mol.bonds.push(Bond::new(1, 3, BondOrder::Single));

        let (charges, h_charges) = gasteiger_charges_with_hydrogens(&mol);
        assert!((charges[2] - charges[3]).abs() < 1e-9);
        let total: f64 = charges.iter().chain(h_charges.iter()).sum();
        assert!((total + 1.0).abs() < 1e-9);
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
