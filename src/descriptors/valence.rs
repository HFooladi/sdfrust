//! Atom valence, degree, and implicit hydrogen count calculations.
//!
//! Provides functions for computing atom degrees and implicit hydrogen counts
//! based on standard valence rules. These are key features for GNN models
//! (2 of 9 OGB atom features: degree, num_hs).
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::valence;
//!
//! let mut mol = Molecule::new("methane");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(3, "H", 0.0, 1.0, 0.0));
//! mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
//!
//! assert_eq!(valence::atom_degree(&mol, 0), 4);
//! assert_eq!(valence::implicit_hydrogen_count(&mol, 0), 0); // all H explicit
//! ```

use crate::descriptors::elements::{element_by_atomic_number, get_element};
use crate::graph::{AdjacencyList, is_hydrogen};
use crate::molecule::Molecule;

/// Allowed valences for an element with the given formal charge, following
/// RDKit: a charged main-group atom uses the valences of the isoelectronic
/// neutral element (atomic number minus charge). Otherwise the neutral list is
/// shifted by the charge (reversed for early elements).
///
/// Returns `None` for unknown elements and `Some(&[])`-like empty results for
/// elements that take any valence (metals), which get no implicit hydrogens.
fn allowed_valences(element: &str, charge: i8) -> Option<Vec<i32>> {
    let data = get_element(element)?;
    if data.atomic_number == 1 && charge != 0 {
        return Some(vec![0]); // H+ / H- carry no further hydrogens
    }
    // A charged main-group atom takes the valences of its isoelectronic
    // neutral element (N+ like C, O- like F, B- like C, B-2 like N, C+ like B)
    if charge != 0 && !data.valences.is_empty() {
        let iso = data.atomic_number as i32 - charge as i32;
        if let Some(iso) = u8::try_from(iso)
            .ok()
            .and_then(element_by_atomic_number)
            .filter(|e| !e.valences.is_empty())
        {
            return Some(iso.valences.iter().map(|&v| v as i32).collect());
        }
    }
    let mut shift = charge as i32;
    if data.outer_electrons < 4 || (data.atomic_number == 6 && charge > 0) {
        shift = -shift;
    }
    Some(
        data.valences
            .iter()
            .map(|&v| v as i32 + shift)
            .filter(|&v| v >= 0)
            .collect(),
    )
}

/// Returns the degree (number of explicit bonds) for atom `idx`.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::valence;
///
/// let mut mol = Molecule::new("ethanol");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
///
/// assert_eq!(valence::atom_degree(&mol, 0), 1);
/// assert_eq!(valence::atom_degree(&mol, 1), 2);
/// ```
pub fn atom_degree(mol: &Molecule, idx: usize) -> usize {
    let adj = AdjacencyList::from_molecule(mol);
    adj.degree(idx)
}

/// Returns the sum of bond orders for atom `idx`.
///
/// Aromatic bonds contribute 1.5 each (rounded up for valence purposes).
pub fn bond_order_sum(mol: &Molecule, idx: usize) -> f64 {
    mol.bonds
        .iter()
        .filter(|b| b.contains_atom(idx))
        .map(|b| b.order.order())
        .sum()
}

/// Returns the integer bond order sum for valence calculation.
///
/// For aromatic bonds, each counts as 1.5 but the total is rounded to nearest integer.
pub fn bond_order_sum_int(mol: &Molecule, idx: usize) -> u8 {
    let sum: f64 = bond_order_sum(mol, idx);
    sum.round() as u8
}

/// Returns the implicit hydrogen count for atom `idx`.
///
/// Follows RDKit's valence model: the result fills the atom up to the smallest
/// allowed valence (for its element and formal charge) that is at least
/// `bond_order_sum + radical_electrons`. Metals and unknown elements get none.
/// The molfile valence field (`vvv`) and hydrogen-count field override this.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::valence;
///
/// // Methanol: C has 1 bond to O, so implicit H = 4 - 1 = 3
/// let mut mol = Molecule::new("methanol_skeleton");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "O", 1.4, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// assert_eq!(valence::implicit_hydrogen_count(&mol, 0), 3); // C with 1 bond → 3 implicit H
/// assert_eq!(valence::implicit_hydrogen_count(&mol, 1), 1); // O with 1 bond → 1 implicit H
/// ```
pub fn implicit_hydrogen_count(mol: &Molecule, idx: usize) -> u8 {
    let atom = match mol.atoms.get(idx) {
        Some(a) => a,
        None => return 0,
    };

    // If the SDF file specified an explicit hydrogen count, use it.
    // SDF hydrogen_count field: 0 = use default, 1 = H0, 2 = H1, etc.
    if let Some(hc) = atom.hydrogen_count {
        if hc > 0 {
            return hc.saturating_sub(1); // SDF convention: 1=H0, 2=H1, ...
        }
    }

    if atom.no_implicit {
        return 0;
    }

    let bonds = bond_order_sum_int(mol, idx) as i32;

    // An explicit valence from the molfile (vvv field) is bonds + hydrogens
    // (radical electrons are not part of it); 15 means zero valence.
    if let Some(v) = atom.valence {
        return if v == 15 {
            0
        } else {
            (v as i32 - bonds).clamp(0, u8::MAX as i32) as u8
        };
    }

    let used = bonds + atom.radical_electrons() as i32;

    // Smallest allowed valence that accommodates the explicit bonds and radicals.
    let Some(valences) = allowed_valences(&atom.element, atom.formal_charge) else {
        return 0;
    };
    valences
        .iter()
        .find(|&&v| v >= used)
        .map(|&v| (v - used) as u8)
        .unwrap_or(0)
}

/// Returns the total hydrogen count (implicit + explicit H neighbors) for atom `idx`.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::valence;
///
/// let mut mol = Molecule::new("methane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
///
/// // C has 2 explicit H + 2 implicit H = 4 total
/// assert_eq!(valence::total_hydrogen_count(&mol, 0), 4);
/// ```
pub fn total_hydrogen_count(mol: &Molecule, idx: usize) -> u8 {
    let implicit = implicit_hydrogen_count(mol, idx);
    let explicit = explicit_hydrogen_count(mol, idx);
    implicit.saturating_add(explicit)
}

/// Returns the number of explicit hydrogen neighbors for atom `idx`.
pub fn explicit_hydrogen_count(mol: &Molecule, idx: usize) -> u8 {
    let adj = AdjacencyList::from_molecule(mol);
    adj.neighbors(idx)
        .iter()
        .filter(|(neighbor, _)| {
            mol.atoms
                .get(*neighbor)
                .map(|a| is_hydrogen(&a.element))
                .unwrap_or(false)
        })
        .count() as u8
}

/// Compute all atom degrees for the molecule.
///
/// Returns a vector of length `mol.atom_count()` with each atom's degree.
pub fn all_atom_degrees(mol: &Molecule) -> Vec<usize> {
    let adj = AdjacencyList::from_molecule(mol);
    adj.degrees().to_vec()
}

/// Compute all implicit hydrogen counts for the molecule.
///
/// Returns a vector of length `mol.atom_count()`.
pub fn all_implicit_hydrogen_counts(mol: &Molecule) -> Vec<u8> {
    (0..mol.atom_count())
        .map(|i| implicit_hydrogen_count(mol, i))
        .collect()
}

/// Compute all total hydrogen counts for the molecule.
///
/// Returns a vector of length `mol.atom_count()`.
pub fn all_total_hydrogen_counts(mol: &Molecule) -> Vec<u8> {
    (0..mol.atom_count())
        .map(|i| total_hydrogen_count(mol, i))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    fn make_ethanol() -> Molecule {
        // C-C-O with explicit hydrogens on one carbon
        let mut mol = Molecule::new("ethanol");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0)); // CH3
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0)); // CH2
        mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0)); // OH
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol
    }

    #[test]
    fn test_atom_degree() {
        let mol = make_ethanol();
        assert_eq!(atom_degree(&mol, 0), 1);
        assert_eq!(atom_degree(&mol, 1), 2);
        assert_eq!(atom_degree(&mol, 2), 1);
    }

    #[test]
    fn test_bond_order_sum_single() {
        let mol = make_ethanol();
        assert!((bond_order_sum(&mol, 0) - 1.0).abs() < 1e-10);
        assert!((bond_order_sum(&mol, 1) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_bond_order_sum_double() {
        let mut mol = Molecule::new("CO");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        assert!((bond_order_sum(&mol, 0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_implicit_h_methanol_skeleton() {
        // C-O with no explicit H
        let mut mol = Molecule::new("methanol_skel");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.4, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        assert_eq!(implicit_hydrogen_count(&mol, 0), 3); // C valence 4 - 1 bond = 3
        assert_eq!(implicit_hydrogen_count(&mol, 1), 1); // O valence 2 - 1 bond = 1
    }

    #[test]
    fn test_implicit_h_formaldehyde() {
        // C=O skeleton
        let mut mol = Molecule::new("formaldehyde_skel");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        assert_eq!(implicit_hydrogen_count(&mol, 0), 2); // C valence 4 - 2 = 2
        assert_eq!(implicit_hydrogen_count(&mol, 1), 0); // O valence 2 - 2 = 0
    }

    #[test]
    fn test_implicit_h_nitrogen() {
        // NH3: just N with no bonds (all H implicit)
        let mut mol = Molecule::new("ammonia");
        mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
        assert_eq!(implicit_hydrogen_count(&mol, 0), 3);
    }

    #[test]
    fn test_implicit_h_charged_nitrogen() {
        // NH4+: N+ with 0 bonds → default valence 4
        let mut mol = Molecule::new("ammonium");
        let mut atom = Atom::new(0, "N", 0.0, 0.0, 0.0);
        atom.formal_charge = 1;
        mol.atoms.push(atom);
        assert_eq!(implicit_hydrogen_count(&mol, 0), 4);
    }

    #[test]
    fn test_total_h_methane_partial() {
        // C with 2 explicit H bonds → implicit 2 + explicit 2 = 4
        let mut mol = Molecule::new("methane_partial");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

        assert_eq!(total_hydrogen_count(&mol, 0), 4);
        assert_eq!(explicit_hydrogen_count(&mol, 0), 2);
        assert_eq!(implicit_hydrogen_count(&mol, 0), 2);
    }

    #[test]
    fn test_implicit_h_aromatic() {
        // Benzene ring: each C has 2 aromatic bonds (sum ~3), so implicit H = 4 - 3 = 1
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }
        // Each C has 2 aromatic bonds = 2 * 1.5 = 3.0, rounded to 3
        // implicit H = 4 - 3 = 1
        assert_eq!(implicit_hydrogen_count(&mol, 0), 1);
    }

    #[test]
    fn test_all_atom_degrees() {
        let mol = make_ethanol();
        let degrees = all_atom_degrees(&mol);
        assert_eq!(degrees, vec![1, 2, 1]);
    }

    #[test]
    fn test_all_implicit_h() {
        let mol = make_ethanol();
        let h_counts = all_implicit_hydrogen_counts(&mol);
        assert_eq!(h_counts, vec![3, 2, 1]); // CH3, CH2, OH
    }

    #[test]
    fn test_default_valence_unknown() {
        assert_eq!(allowed_valences("Xx", 0), None);
    }

    #[test]
    fn test_sdf_hydrogen_count_field() {
        let mut mol = Molecule::new("test");
        let mut atom = Atom::new(0, "C", 0.0, 0.0, 0.0);
        atom.hydrogen_count = Some(2); // SDF convention: 2 = H1
        mol.atoms.push(atom);
        // Should use the explicit SDF field (2 - 1 = 1)
        assert_eq!(implicit_hydrogen_count(&mol, 0), 1);
    }
}
