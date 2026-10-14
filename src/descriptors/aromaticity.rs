//! Aromaticity perception for molecular rings.
//!
//! Follows RDKit's default model:
//! 1. Bonds annotated as aromatic in the file (`BondOrder::Aromatic`) are trusted.
//! 2. Every ring atom gets a π-electron count from its own bonds and charge
//!    (an atom double-bonded to another ring atom gives 1, pyrrole-type N and
//!    furan/thiophene O/S give 2, an exocyclic C=O leaves the carbon with 0, ...).
//! 3. Hückel's 4n+2 rule is tested on each SSSR ring and on every union of
//!    fused rings, so Kekulé naphthalene, indole and azulene are recognized.
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
use crate::descriptors::elements::get_element;
use crate::descriptors::rings::{Ring, sssr};
use crate::descriptors::valence::implicit_hydrogen_count;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// Fused ring systems up to this many rings have every union tested;
/// larger systems test unions of at most three rings.
const MAX_EXHAUSTIVE_RINGS: usize = 10;

/// Pauling electronegativity for the elements that matter for exocyclic
/// double bonds; unknown elements count as 2.0.
fn electronegativity(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 2.20,
        5 => 2.04,
        6 => 2.55,
        7 => 3.04,
        8 => 3.44,
        9 => 3.98,
        14 => 1.90,
        15 => 2.19,
        16 => 2.58,
        17 => 3.16,
        33 => 2.18,
        34 => 2.55,
        35 => 2.96,
        52 => 2.10,
        53 => 2.66,
        _ => 2.0,
    }
}

/// π electrons that ring atom `idx` can donate to an aromatic system, or
/// `None` if the atom cannot be part of one (e.g. an sp3 carbon).
fn pi_electrons(mol: &Molecule, adj: &AdjacencyList, ring_bond: &[bool], idx: usize) -> Option<u8> {
    let atom = &mol.atoms[idx];
    let z = get_element(&atom.element)?.atomic_number;
    let total_degree = adj.degree(idx) + implicit_hydrogen_count(mol, idx) as usize;
    if total_degree > 3 {
        return None;
    }

    let mut ring_pi = false;
    let mut exocyclic_partner = None;
    for &(neighbor, bond_idx) in adj.neighbors(idx) {
        match mol.bonds[bond_idx].order {
            BondOrder::Aromatic => ring_pi = true,
            BondOrder::Double if ring_bond[bond_idx] => ring_pi = true,
            BondOrder::Double => exocyclic_partner = Some(neighbor),
            BondOrder::Triple => return None,
            _ => {}
        }
    }
    if ring_pi {
        return Some(1);
    }
    if let Some(partner) = exocyclic_partner {
        // A ring carbon with an exocyclic double bond to a more electronegative
        // atom (C=O, C=N, C=S) contributes an empty p orbital; any other
        // exocyclic double bond (including S=O, as in sulfoxides) breaks the ring.
        let partner_z = get_element(&mol.atoms[partner].element)?.atomic_number;
        return (z == 6 && electronegativity(partner_z) > electronegativity(z)).then_some(0);
    }

    // No double bond: lone pairs, empty orbitals or charges decide.
    match (z, atom.formal_charge) {
        (6, -1) => Some(2),                                    // cyclopentadienide
        (6, 1) => Some(0),                                     // tropylium
        (7 | 15 | 33, 0) if total_degree == 3 => Some(2),      // pyrrole-type
        (7 | 15 | 33, -1) if total_degree == 2 => Some(2),     // pyrrolide
        (8 | 16 | 34 | 52, 0) if total_degree == 2 => Some(2), // furan, thiophene
        (5, 0) if total_degree == 3 => Some(0),                // borole-type vacant
        _ => None,
    }
}

/// Aromatic flags for every atom and bond.
struct Aromaticity {
    atoms: Vec<bool>,
    bonds: Vec<bool>,
}

/// Runs aromaticity perception once for the whole molecule.
fn perceive(mol: &Molecule) -> Aromaticity {
    let n = mol.atom_count();
    let m = mol.bond_count();
    let mut result = Aromaticity {
        atoms: vec![false; n],
        bonds: vec![false; m],
    };

    // Stage 1: file-annotated aromatic bonds
    for (i, bond) in mol.bonds.iter().enumerate() {
        if bond.order == BondOrder::Aromatic && bond.atom1 < n && bond.atom2 < n {
            result.bonds[i] = true;
            result.atoms[bond.atom1] = true;
            result.atoms[bond.atom2] = true;
        }
    }

    // Stage 2: Hückel perception on rings and fused ring unions
    let rings = sssr(mol);
    if rings.is_empty() {
        return result;
    }
    let adj = AdjacencyList::from_molecule(mol);
    let mut ring_bond = vec![false; m];
    for ring in &rings {
        for &b in &ring.bonds {
            ring_bond[b] = true;
        }
    }
    let mut electrons: Vec<Option<Option<u8>>> = vec![None; n];
    let mut donation = |idx: usize| -> Option<u8> {
        *electrons[idx].get_or_insert_with(|| pi_electrons(mol, &adj, &ring_bond, idx))
    };

    // Only rings whose atoms can all donate electrons can be (part of) aromatic systems
    let candidates: Vec<&Ring> = rings
        .iter()
        .filter(|r| r.atoms.iter().all(|&a| donation(a).is_some()))
        .collect();
    let electrons: Vec<u8> = (0..n)
        .map(|i| electrons[i].flatten().unwrap_or(0))
        .collect();

    for system in fused_systems(&candidates) {
        for subset in ring_unions(&candidates, &system) {
            let mut atom_in = vec![false; n];
            for &r in &subset {
                for &a in &candidates[r].atoms {
                    atom_in[a] = true;
                }
            }
            let total: u32 = (0..n)
                .filter(|&a| atom_in[a])
                .map(|a| electrons[a] as u32)
                .sum();
            if total >= 2 && (total - 2) % 4 == 0 {
                // All atoms are aromatic, but for a union of fused rings only the
                // perimeter bonds are (a fusion bond needs an aromatic ring of its own).
                let mut bond_uses = std::collections::HashMap::new();
                for &r in &subset {
                    for &a in &candidates[r].atoms {
                        result.atoms[a] = true;
                    }
                    for &b in &candidates[r].bonds {
                        *bond_uses.entry(b).or_insert(0) += 1;
                    }
                }
                for (b, uses) in bond_uses {
                    if uses == 1 {
                        result.bonds[b] = true;
                    }
                }
            }
        }
    }

    result
}

/// Groups rings into fused systems (rings connected through shared bonds).
fn fused_systems(rings: &[&Ring]) -> Vec<Vec<usize>> {
    let shares_bond = |a: &Ring, b: &Ring| a.bonds.iter().any(|x| b.bonds.contains(x));
    let mut system_of: Vec<Option<usize>> = vec![None; rings.len()];
    let mut systems: Vec<Vec<usize>> = Vec::new();
    for start in 0..rings.len() {
        if system_of[start].is_some() {
            continue;
        }
        let id = systems.len();
        let mut members = vec![start];
        system_of[start] = Some(id);
        let mut i = 0;
        while i < members.len() {
            let current = members[i];
            for other in 0..rings.len() {
                if system_of[other].is_none() && shares_bond(rings[current], rings[other]) {
                    system_of[other] = Some(id);
                    members.push(other);
                }
            }
            i += 1;
        }
        systems.push(members);
    }
    systems
}

/// Connected subsets of a fused system's rings to test for aromaticity.
fn ring_unions(rings: &[&Ring], system: &[usize]) -> Vec<Vec<usize>> {
    let k = system.len();
    let max_size = if k <= MAX_EXHAUSTIVE_RINGS { k } else { 3 };
    let adjacent = |a: usize, b: usize| rings[a].bonds.iter().any(|x| rings[b].bonds.contains(x));
    let mut unions = Vec::new();
    let mut frontier: Vec<Vec<usize>> = system.iter().map(|&r| vec![r]).collect();
    let mut seen = std::collections::HashSet::new();
    for size in 1..=max_size {
        let mut next = Vec::new();
        for subset in frontier {
            let mut key = subset.clone();
            key.sort_unstable();
            if !seen.insert(key.clone()) {
                continue;
            }
            if size < max_size {
                for &r in system {
                    if !subset.contains(&r) && subset.iter().any(|&s| adjacent(s, r)) {
                        let mut grown = subset.clone();
                        grown.push(r);
                        next.push(grown);
                    }
                }
            }
            unions.push(key);
        }
        frontier = next;
    }
    unions
}

/// Check if an atom is aromatic.
///
/// An atom is aromatic if it is bonded by an aromatic bond in the file, or is
/// part of a ring (or fused ring union) satisfying Hückel's 4n+2 rule. For many
/// atoms, prefer [`all_aromatic_atoms`], which perceives the molecule once.
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
    idx < mol.atoms.len() && perceive(mol).atoms[idx]
}

/// Check if a bond is aromatic (annotated in the file or part of an aromatic
/// ring system).
pub fn is_aromatic_bond(mol: &Molecule, bond_idx: usize) -> bool {
    bond_idx < mol.bonds.len() && perceive(mol).bonds[bond_idx]
}

/// Compute aromaticity for all atoms.
///
/// Returns a vector of booleans of length `mol.atom_count()`.
pub fn all_aromatic_atoms(mol: &Molecule) -> Vec<bool> {
    perceive(mol).atoms
}

/// Compute aromaticity for all bonds.
///
/// Returns a vector of booleans of length `mol.bond_count()`.
pub fn all_aromatic_bonds(mol: &Molecule) -> Vec<bool> {
    perceive(mol).bonds
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
