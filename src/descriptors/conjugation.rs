//! Conjugation detection for bonds.
//!
//! A bond is conjugated if:
//! - It is an aromatic bond, OR
//! - It is a single bond between two SP2 atoms, OR
//! - It is a double bond adjacent to another double or aromatic bond
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
//! assert!(conj[1]); // C-C between two SP2 atoms is conjugated
//! assert!(conj[2]); // C=C is conjugated
//! ```

use crate::bond::BondOrder;
use crate::descriptors::hybridization::{Hybridization, all_hybridizations};
use crate::molecule::Molecule;

/// Check if a bond is conjugated.
///
/// A bond is conjugated if any of these conditions hold:
/// - It is aromatic
/// - It is a single bond and both atoms are SP2 hybridized
/// - It is a double bond adjacent to another double or aromatic bond
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
    let bond = match mol.bonds.get(bond_idx) {
        Some(b) => b,
        None => return false,
    };

    // Rule 1: Aromatic bonds are always conjugated
    if bond.is_aromatic() {
        return true;
    }

    let hybs = all_hybridizations(mol);

    // Rule 2: Single bond between two SP2 atoms
    if bond.order == BondOrder::Single {
        let hyb1 = hybs
            .get(bond.atom1)
            .copied()
            .unwrap_or(Hybridization::Other);
        let hyb2 = hybs
            .get(bond.atom2)
            .copied()
            .unwrap_or(Hybridization::Other);
        if hyb1 == Hybridization::SP2 && hyb2 == Hybridization::SP2 {
            return true;
        }
    }

    // Rule 3: Double/triple bond adjacent to unsaturation
    if bond.order == BondOrder::Double || bond.order == BondOrder::Triple {
        let hybs = all_hybridizations(mol);
        let hyb1 = hybs
            .get(bond.atom1)
            .copied()
            .unwrap_or(Hybridization::Other);
        let hyb2 = hybs
            .get(bond.atom2)
            .copied()
            .unwrap_or(Hybridization::Other);
        let is_sp2_or_sp = |h: Hybridization| h == Hybridization::SP2 || h == Hybridization::SP;

        if is_sp2_or_sp(hyb1) && is_sp2_or_sp(hyb2) {
            for (other_idx, other_bond) in mol.bonds.iter().enumerate() {
                if other_idx == bond_idx {
                    continue;
                }
                let shares_atom =
                    other_bond.contains_atom(bond.atom1) || other_bond.contains_atom(bond.atom2);
                if shares_atom
                    && (other_bond.order == BondOrder::Double
                        || other_bond.order == BondOrder::Triple
                        || other_bond.is_aromatic()
                        || {
                            if other_bond.order == BondOrder::Single {
                                let other_end = if other_bond.contains_atom(bond.atom1) {
                                    other_bond.other_atom(bond.atom1).unwrap_or(usize::MAX)
                                } else {
                                    other_bond.other_atom(bond.atom2).unwrap_or(usize::MAX)
                                };
                                let other_hyb =
                                    hybs.get(other_end).copied().unwrap_or(Hybridization::Other);
                                is_sp2_or_sp(other_hyb)
                            } else {
                                false
                            }
                        })
                {
                    return true;
                }
            }
        }
    }

    false
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
    let m = mol.bond_count();
    if m == 0 {
        return vec![];
    }

    let hybs = all_hybridizations(mol);
    let mut result = vec![false; m];

    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        // Rule 1: Aromatic
        if bond.is_aromatic() {
            result[bond_idx] = true;
            continue;
        }

        // Rule 2: Single bond between SP2 atoms
        if bond.order == BondOrder::Single {
            let hyb1 = hybs
                .get(bond.atom1)
                .copied()
                .unwrap_or(Hybridization::Other);
            let hyb2 = hybs
                .get(bond.atom2)
                .copied()
                .unwrap_or(Hybridization::Other);
            if hyb1 == Hybridization::SP2 && hyb2 == Hybridization::SP2 {
                result[bond_idx] = true;
                continue;
            }
        }

        // Rule 3: Double/triple bond where both atoms are SP2 (or SP)
        // This covers cases like butadiene C=C-C=C where each C=C bond
        // has SP2 atoms on both sides.
        if bond.order == BondOrder::Double || bond.order == BondOrder::Triple {
            let hyb1 = hybs
                .get(bond.atom1)
                .copied()
                .unwrap_or(Hybridization::Other);
            let hyb2 = hybs
                .get(bond.atom2)
                .copied()
                .unwrap_or(Hybridization::Other);
            let is_sp2_or_sp = |h: Hybridization| h == Hybridization::SP2 || h == Hybridization::SP;
            if is_sp2_or_sp(hyb1) && is_sp2_or_sp(hyb2) {
                // Check if either atom has another unsaturated neighbor
                let mut is_conj = false;
                for (other_idx, other_bond) in mol.bonds.iter().enumerate() {
                    if other_idx == bond_idx {
                        continue;
                    }
                    let shares_atom = other_bond.contains_atom(bond.atom1)
                        || other_bond.contains_atom(bond.atom2);
                    if shares_atom
                        && (other_bond.order == BondOrder::Double
                            || other_bond.order == BondOrder::Triple
                            || other_bond.is_aromatic()
                            || {
                                // Single bond to another SP2 atom (conjugation chain)
                                if other_bond.order == BondOrder::Single {
                                    let other_end = if other_bond.contains_atom(bond.atom1) {
                                        other_bond.other_atom(bond.atom1).unwrap_or(usize::MAX)
                                    } else {
                                        other_bond.other_atom(bond.atom2).unwrap_or(usize::MAX)
                                    };
                                    let other_hyb = hybs
                                        .get(other_end)
                                        .copied()
                                        .unwrap_or(Hybridization::Other);
                                    is_sp2_or_sp(other_hyb)
                                } else {
                                    false
                                }
                            })
                    {
                        is_conj = true;
                        break;
                    }
                }
                if is_conj {
                    result[bond_idx] = true;
                }
            }
        }
    }

    result
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
