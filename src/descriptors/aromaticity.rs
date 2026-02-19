//! Aromaticity perception for molecular rings.
//!
//! Two-stage aromaticity:
//! 1. Trust existing `BondOrder::Aromatic` from the parsed file
//! 2. Perceive aromaticity from Kekulized form via Hückel 4n+2 rule on SSSR rings
//!
//! Provides 1 of 9 OGB atom features (is_aromatic).
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::aromaticity;
//!
//! let mut mol = Molecule::new("benzene");
//! for i in 0..6 {
//!     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
//! }
//! for i in 0..6 {
//!     mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
//! }
//!
//! assert!(aromaticity::is_aromatic_atom(&mol, 0));
//! ```

use crate::bond::BondOrder;
use crate::descriptors::rings::sssr;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Count pi electrons contributed by an atom in a ring.
///
/// Uses a lookup table for common heteroatoms in aromatic rings.
///
/// # Arguments
/// * `in_double_bond` - atom has a double bond to a ring neighbor
/// * `has_exocyclic_double` - atom has a double bond to a non-ring neighbor (e.g., C=O)
fn pi_electrons(
    element: &str,
    charge: i8,
    bond_order_sum: f64,
    in_double_bond: bool,
    has_exocyclic_double: bool,
) -> Option<u8> {
    let elem = element.trim();
    // Normalize: uppercase first, lowercase rest
    let upper: String = elem
        .chars()
        .next()
        .map(|c| c.to_uppercase().collect::<String>())
        .unwrap_or_default()
        + &elem.chars().skip(1).collect::<String>().to_lowercase();

    match upper.as_str() {
        "C" => {
            if in_double_bond {
                Some(1) // C in ring C=C contributes 1 pi electron
            } else if charge == -1 {
                Some(2) // C- (carbanion in ring) can contribute 2
            } else if has_exocyclic_double {
                Some(0) // sp2 C with exocyclic double bond (e.g., C=O in purines)
            } else {
                None // sp3 C with no double bonds breaks conjugation
            }
        }
        "N" => {
            if in_double_bond {
                Some(1) // Pyridine-like N (=N-)
            } else if charge == 0 && bond_order_sum <= 3.0 {
                Some(2) // Pyrrole-like N (-NH-)
            } else if charge == 1 {
                Some(1) // N+ in ring
            } else {
                Some(0)
            }
        }
        "O" => {
            if in_double_bond {
                Some(1)
            } else {
                Some(2) // Furan-like O
            }
        }
        "S" => {
            if in_double_bond {
                Some(1)
            } else {
                Some(2) // Thiophene-like S
            }
        }
        "Se" => {
            if in_double_bond {
                Some(1)
            } else {
                Some(2)
            }
        }
        "P" => {
            if in_double_bond {
                Some(1)
            } else {
                Some(2)
            }
        }
        _ => None,
    }
}

/// Check if a ring satisfies the Hückel 4n+2 rule.
///
/// A ring is aromatic if:
/// - The ring contains at least one double or aromatic bond (all-single-bond rings are never aromatic)
/// - All atoms can contribute a defined number of pi electrons
/// - The total pi electron count satisfies 4n+2 (n=0,1,2,...)
fn is_huckel_aromatic(ring_atoms: &[usize], mol: &Molecule, adj: &AdjacencyList) -> bool {
    let ring_set: std::collections::HashSet<usize> = ring_atoms.iter().copied().collect();

    // First pass: check that the ring has at least one double/aromatic bond.
    // Fully saturated rings (like glucose pyranose) should never be aromatic.
    let mut has_pi_bond = false;
    for &atom_idx in ring_atoms {
        for &(neighbor, bond_idx) in adj.neighbors(atom_idx) {
            if ring_set.contains(&neighbor) {
                if let Some(bond) = mol.bonds.get(bond_idx) {
                    if matches!(
                        bond.order,
                        BondOrder::Double
                            | BondOrder::Triple
                            | BondOrder::Aromatic
                            | BondOrder::DoubleOrAromatic
                    ) {
                        has_pi_bond = true;
                        break;
                    }
                }
            }
        }
        if has_pi_bond {
            break;
        }
    }
    if !has_pi_bond {
        return false;
    }

    // Second pass: count pi electrons
    let mut total_pi = 0u16;
    for &atom_idx in ring_atoms {
        let atom = match mol.atoms.get(atom_idx) {
            Some(a) => a,
            None => return false,
        };

        // Check if this atom has a double bond to another ring atom,
        // and/or an exocyclic double bond (to a non-ring atom)
        let mut in_ring_double = false;
        let mut has_exocyclic_double = false;
        let bo_sum: f64 = adj
            .neighbors(atom_idx)
            .iter()
            .filter_map(|&(neighbor, bond_idx)| {
                let bond = mol.bonds.get(bond_idx)?;
                if matches!(
                    bond.order,
                    BondOrder::Double | BondOrder::Aromatic | BondOrder::DoubleOrAromatic
                ) {
                    if ring_set.contains(&neighbor) {
                        in_ring_double = true;
                    } else {
                        has_exocyclic_double = true;
                    }
                }
                Some(bond.order.order())
            })
            .sum();

        match pi_electrons(
            &atom.element,
            atom.formal_charge,
            bo_sum,
            in_ring_double,
            has_exocyclic_double,
        ) {
            Some(pe) => total_pi += pe as u16,
            None => return false, // Unknown atom type or sp3 carbon breaks ring
        }
    }

    // Hückel rule: 4n+2 for n=0,1,2,...
    // Valid values: 2, 6, 10, 14, 18, 22, ...
    if total_pi < 2 {
        return false;
    }
    (total_pi - 2) % 4 == 0
}

/// Check if an atom is aromatic.
///
/// An atom is aromatic if:
/// - It is bonded by at least one aromatic bond (from file), OR
/// - It is part of a ring that satisfies the Hückel 4n+2 rule
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::aromaticity;
///
/// let mut mol = Molecule::new("benzene");
/// for i in 0..6 {
///     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
///
/// for i in 0..6 {
///     assert!(aromaticity::is_aromatic_atom(&mol, i));
/// }
/// ```
pub fn is_aromatic_atom(mol: &Molecule, idx: usize) -> bool {
    if idx >= mol.atoms.len() {
        return false;
    }

    // Stage 1: Trust file-annotated aromatic bonds
    for bond in &mol.bonds {
        if bond.contains_atom(idx) && bond.is_aromatic() {
            return true;
        }
    }

    // Stage 2: Perceive from Kekulized form using SSSR + Hückel
    let rings = sssr(mol);
    let adj = AdjacencyList::from_molecule(mol);

    for ring in &rings {
        if ring.contains_atom(idx) && is_huckel_aromatic(&ring.atoms, mol, &adj) {
            return true;
        }
    }

    false
}

/// Check if a bond is aromatic.
///
/// A bond is aromatic if:
/// - Its order is `BondOrder::Aromatic` (from file), OR
/// - Both its atoms are in the same aromatic ring
pub fn is_aromatic_bond(mol: &Molecule, bond_idx: usize) -> bool {
    let bond = match mol.bonds.get(bond_idx) {
        Some(b) => b,
        None => return false,
    };

    // Stage 1: trust file
    if bond.is_aromatic() {
        return true;
    }

    // Stage 2: both atoms in same aromatic ring
    let rings = sssr(mol);
    let adj = AdjacencyList::from_molecule(mol);

    for ring in &rings {
        if ring.contains_atom(bond.atom1)
            && ring.contains_atom(bond.atom2)
            && is_huckel_aromatic(&ring.atoms, mol, &adj)
        {
            return true;
        }
    }

    false
}

/// Compute aromaticity for all atoms.
///
/// Returns a vector of booleans of length `mol.atom_count()`.
pub fn all_aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut result = vec![false; n];

    // Stage 1: file-annotated aromatic bonds
    for bond in &mol.bonds {
        if bond.is_aromatic() {
            if bond.atom1 < n {
                result[bond.atom1] = true;
            }
            if bond.atom2 < n {
                result[bond.atom2] = true;
            }
        }
    }

    // Stage 2: Hückel perception for atoms not yet marked
    if result.iter().all(|&x| x) {
        return result;
    }

    let rings = sssr(mol);
    let adj = AdjacencyList::from_molecule(mol);

    for ring in &rings {
        if is_huckel_aromatic(&ring.atoms, mol, &adj) {
            for &atom_idx in &ring.atoms {
                if atom_idx < n {
                    result[atom_idx] = true;
                }
            }
        }
    }

    result
}

/// Compute aromaticity for all bonds.
///
/// Returns a vector of booleans of length `mol.bond_count()`.
pub fn all_aromatic_bonds(mol: &Molecule) -> Vec<bool> {
    let m = mol.bond_count();
    let mut result = vec![false; m];

    // Stage 1: file-annotated
    for (i, bond) in mol.bonds.iter().enumerate() {
        if bond.is_aromatic() {
            result[i] = true;
        }
    }

    // Stage 2: Hückel perception
    let rings = sssr(mol);
    let adj = AdjacencyList::from_molecule(mol);

    for ring in &rings {
        if is_huckel_aromatic(&ring.atoms, mol, &adj) {
            for &bond_idx in &ring.bonds {
                if bond_idx < m {
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

    fn make_benzene_aromatic() -> Molecule {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }
        mol
    }

    fn make_benzene_kekulized() -> Molecule {
        // Alternating single/double bonds (Kekulé form)
        let mut mol = Molecule::new("benzene_kekule");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            mol.bonds.push(Bond::new(i, (i + 1) % 6, order));
        }
        mol
    }

    #[test]
    fn test_benzene_aromatic_bonds() {
        let mol = make_benzene_aromatic();
        for i in 0..6 {
            assert!(is_aromatic_atom(&mol, i));
        }
    }

    #[test]
    fn test_benzene_kekulized_perception() {
        let mol = make_benzene_kekulized();
        // Should perceive aromaticity via Hückel rule
        // 6 C atoms each contributing 1 pi electron = 6 total → 4n+2 for n=1
        for i in 0..6 {
            assert!(
                is_aromatic_atom(&mol, i),
                "Atom {} should be aromatic in Kekulized benzene",
                i
            );
        }
    }

    #[test]
    fn test_cyclopentane_not_aromatic() {
        // 5 sp3 carbons with single bonds → not aromatic
        let mut mol = Molecule::new("cyclopentane");
        for i in 0..5 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..5 {
            mol.bonds.push(Bond::new(i, (i + 1) % 5, BondOrder::Single));
        }

        for i in 0..5 {
            assert!(!is_aromatic_atom(&mol, i));
        }
    }

    #[test]
    fn test_pyrrole() {
        // Pyrrole: 4C + 1N in 5-membered ring with alternating bonds
        // N contributes 2 pi electrons, each C=C pair contributes 2 → total 6 → 4n+2
        let mut mol = Molecule::new("pyrrole");
        mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 1.5, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 0.5, 1.5, 0.0));
        mol.atoms.push(Atom::new(4, "C", -0.5, 1.0, 0.0));

        // N-C single, C=C, C-C, C=C, C-N
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(3, 4, BondOrder::Double));
        mol.bonds.push(Bond::new(4, 0, BondOrder::Single));

        // Pyrrole is aromatic (6 pi electrons)
        assert!(is_aromatic_atom(&mol, 0), "N should be aromatic in pyrrole");
    }

    #[test]
    fn test_furan() {
        // Furan: 4C + 1O in 5-membered ring
        let mut mol = Molecule::new("furan");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 1.5, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 0.5, 1.5, 0.0));
        mol.atoms.push(Atom::new(4, "C", -0.5, 1.0, 0.0));

        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(3, 4, BondOrder::Double));
        mol.bonds.push(Bond::new(4, 0, BondOrder::Single));

        assert!(is_aromatic_atom(&mol, 0), "O should be aromatic in furan");
    }

    #[test]
    fn test_ethane_not_aromatic() {
        let mut mol = Molecule::new("ethane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        assert!(!is_aromatic_atom(&mol, 0));
    }

    #[test]
    fn test_all_aromatic_atoms_benzene() {
        let mol = make_benzene_aromatic();
        let arom = all_aromatic_atoms(&mol);
        assert!(arom.iter().all(|&x| x));
    }

    #[test]
    fn test_all_aromatic_bonds_benzene() {
        let mol = make_benzene_aromatic();
        let arom = all_aromatic_bonds(&mol);
        assert!(arom.iter().all(|&x| x));
    }

    #[test]
    fn test_is_aromatic_bond() {
        let mol = make_benzene_aromatic();
        for i in 0..6 {
            assert!(is_aromatic_bond(&mol, i));
        }
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new("empty");
        assert!(!is_aromatic_atom(&mol, 0));
        let arom = all_aromatic_atoms(&mol);
        assert!(arom.is_empty());
    }

    #[test]
    fn test_saturated_ring_not_aromatic() {
        // Tetrahydropyran: 5C + 1O ring, all single bonds (like glucose pyranose)
        let mut mol = Molecule::new("thp");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 1.5, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 0.5, 1.5, 0.0));
        mol.atoms.push(Atom::new(4, "C", -0.5, 1.0, 0.0));
        mol.atoms.push(Atom::new(5, "C", -1.0, 0.0, 0.0));

        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(3, 4, BondOrder::Single));
        mol.bonds.push(Bond::new(4, 5, BondOrder::Single));
        mol.bonds.push(Bond::new(5, 0, BondOrder::Single));

        for i in 0..6 {
            assert!(
                !is_aromatic_atom(&mol, i),
                "Atom {} in saturated ring should not be aromatic",
                i
            );
        }
    }
}
