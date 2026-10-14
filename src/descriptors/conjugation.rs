//! Conjugation detection for bonds.
//!
//! Follows RDKit's rules (`MolOps::setConjugation`):
//! - Aromatic bonds are conjugated.
//! - For an atom with 2 or 3 connections (including hydrogens) that carries a
//!   double, triple or aromatic bond, that bond and each of its other bonds are
//!   conjugated when the neighbour on the other side has at most 3
//!   connections. Both atoms must be able to take part in a π system (have a
//!   π bond or a lone pair); atoms from groups 15/16 below the first row
//!   (P, S, ...) only qualify when terminal, as in RDKit.
//!
//! Completes all 3 OGB bond features (bond_type, stereo, is_conjugated).
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::conjugation;
//!
//! // Butadiene: C=C-C=C (the single bond in the middle is conjugated)
//! let mut mol = Molecule::new("butadiene");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "C", 2.5, 0.0, 0.0));
//! mol.atoms.push(Atom::new(3, "C", 3.8, 0.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
//! mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
//! mol.bonds.push(Bond::new(2, 3, BondOrder::Double));
//!
//! let conj = conjugation::all_conjugated_bonds(&mol);
//! assert!(conj[0]); // C=C is conjugated
//! assert!(conj[1]); // C-C between two C=C is conjugated
//! assert!(conj[2]); // C=C is conjugated
//! ```

use crate::bond::BondOrder;
use crate::descriptors::aromaticity::all_aromatic_bonds;
use crate::descriptors::elements::get_element;
use crate::descriptors::valence::all_implicit_hydrogen_counts;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Per-atom quantities shared by conjugation and hybridization perception.
pub(crate) struct AtomEnvironment {
    pub adj: AdjacencyList,
    /// Implicit hydrogen count per atom.
    pub implicit_h: Vec<u8>,
    /// Sum of explicit bond orders per atom (aromatic = 1.5, rounded).
    pub bond_valence: Vec<i32>,
    /// Perceived aromaticity per bond.
    pub aromatic_bonds: Vec<bool>,
}

impl AtomEnvironment {
    pub fn new(mol: &Molecule) -> Self {
        let adj = AdjacencyList::from_molecule(mol);
        let bond_valence = (0..mol.atom_count())
            .map(|i| {
                adj.neighbors(i)
                    .iter()
                    .map(|&(_, b)| mol.bonds[b].order.order())
                    .sum::<f64>()
                    .round() as i32
            })
            .collect();
        Self {
            implicit_h: all_implicit_hydrogen_counts(mol),
            bond_valence,
            aromatic_bonds: all_aromatic_bonds(mol),
            adj,
        }
    }

    /// Explicit neighbours plus implicit hydrogens.
    pub fn total_degree(&self, idx: usize) -> i32 {
        self.adj.degree(idx) as i32 + self.implicit_h[idx] as i32
    }

    /// Bond orders plus implicit hydrogens.
    pub fn total_valence(&self, idx: usize) -> i32 {
        self.bond_valence[idx] + self.implicit_h[idx] as i32
    }
}

/// Electrons atom `idx` can donate to a π system (RDKit's `countAtomElec`);
/// negative if it cannot take part in one.
fn count_atom_electrons(mol: &Molecule, env: &AtomEnvironment, idx: usize) -> i32 {
    let atom = &mol.atoms[idx];
    let Some(data) = get_element(&atom.element) else {
        return -1;
    };
    let default_valence = data.valences.first().map_or(-1, |&v| v as i32);
    if default_valence <= 1 {
        return -1; // univalent elements can't be conjugated
    }
    let degree = env.total_degree(idx);
    if degree > 3 {
        return -1;
    }
    let lone_pair_electrons =
        (data.outer_electrons as i32 - default_valence - atom.formal_charge as i32).max(0);
    let mut electrons =
        (default_valence - degree) + lone_pair_electrons - atom.radical_electrons() as i32;
    let unsaturations = env.bond_valence[idx] - env.adj.degree(idx) as i32;
    if electrons > 1 && unsaturations > 1 {
        electrons = 1; // a triple bond only offers one electron to the system
    }
    electrons
}

/// Whether atom `idx` can extend a conjugated system (RDKit's `isAtomConjugCand`).
fn is_conjugation_candidate(mol: &Molecule, env: &AtomEnvironment, idx: usize) -> bool {
    let Some(data) = get_element(&mol.atoms[idx].element) else {
        return false;
    };
    let outer = data.outer_electrons;
    let element_ok = data.atomic_number <= 10
        || (outer != 5 && outer != 6)
        || (outer == 6 && env.total_degree(idx) < 2);
    element_ok && count_atom_electrons(mol, env, idx) > 0
}

/// Conjugation flags for every bond, using precomputed atom environments.
pub(crate) fn conjugated_bonds(mol: &Molecule, env: &AtomEnvironment) -> Vec<bool> {
    let mut conjugated = env.aromatic_bonds.clone();
    let is_unsaturated = |b: usize| {
        env.aromatic_bonds[b]
            || matches!(
                mol.bonds[b].order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
            )
    };

    for atom in 0..mol.atom_count() {
        let degree = env.total_degree(atom);
        if !(2..=3).contains(&degree) || !is_conjugation_candidate(mol, env, atom) {
            continue;
        }
        let bonds = env.adj.neighbors(atom);
        for &(_, b1) in bonds {
            if !is_unsaturated(b1) {
                continue;
            }
            for &(other, b2) in bonds {
                if b1 == b2 || env.total_degree(other) > 3 {
                    continue;
                }
                if is_conjugation_candidate(mol, env, other) {
                    conjugated[b1] = true;
                    conjugated[b2] = true;
                }
            }
        }
    }
    conjugated
}

/// Check if a bond is conjugated.
///
/// For many bonds, prefer [`all_conjugated_bonds`], which perceives the
/// molecule once.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::conjugation;
///
/// let mut mol = Molecule::new("benzene");
/// for i in 0..6 {
///     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
///
/// assert!(conjugation::is_conjugated_bond(&mol, 0));
/// ```
pub fn is_conjugated_bond(mol: &Molecule, bond_idx: usize) -> bool {
    bond_idx < mol.bond_count() && all_conjugated_bonds(mol)[bond_idx]
}

/// Compute conjugation for all bonds in the molecule.
///
/// Returns a vector of booleans of length `mol.bond_count()`.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::conjugation;
///
/// let mut mol = Molecule::new("ethane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let conj = conjugation::all_conjugated_bonds(&mol);
/// assert!(!conj[0]); // C-C single bond between SP3 atoms
/// ```
pub fn all_conjugated_bonds(mol: &Molecule) -> Vec<bool> {
    if mol.bond_count() == 0 {
        return vec![];
    }
    conjugated_bonds(mol, &AtomEnvironment::new(mol))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    #[test]
    fn test_benzene_conjugated() {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }

        let conj = all_conjugated_bonds(&mol);
        assert!(conj.iter().all(|&c| c));
    }

    #[test]
    fn test_butadiene() {
        // C=C-C=C: all bonds are conjugated
        let mut mol = Molecule::new("butadiene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 2.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 3.8, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Double));

        let conj = all_conjugated_bonds(&mol);
        assert!(conj[0], "First C=C should be conjugated");
        assert!(conj[1], "C-C between SP2 atoms should be conjugated");
        assert!(conj[2], "Second C=C should be conjugated");
    }

    #[test]
    fn test_ethane_not_conjugated() {
        let mut mol = Molecule::new("ethane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let conj = all_conjugated_bonds(&mol);
        assert!(!conj[0]);
    }

    #[test]
    fn test_isolated_double_bond() {
        // C=C with no adjacent double bonds
        let mut mol = Molecule::new("ethylene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        let conj = all_conjugated_bonds(&mol);
        assert!(!conj[0], "Isolated C=C should not be conjugated");
    }

    #[test]
    fn test_empty() {
        let mol = Molecule::new("empty");
        let conj = all_conjugated_bonds(&mol);
        assert!(conj.is_empty());
    }

    #[test]
    fn test_is_conjugated_bond() {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }

        assert!(is_conjugated_bond(&mol, 0));
    }
}
