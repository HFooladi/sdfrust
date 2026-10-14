//! Smallest Set of Smallest Rings (SSSR) perception.
//!
//! Implements ring perception for molecular graphs, finding the minimum
//! set of independent rings. This is a prerequisite for aromaticity
//! detection and provides ring-based features for ML models.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::descriptors::rings;
//!
//! // Benzene: one 6-membered ring
//! let mut mol = Molecule::new("benzene");
//! for i in 0..6 {
//!     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
//! }
//! for i in 0..6 {
//!     mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
//! }
//!
//! let ring_set = rings::sssr(&mol);
//! assert_eq!(ring_set.len(), 1);
//! assert_eq!(ring_set[0].size(), 6);
//! ```

use crate::graph::AdjacencyList;
use crate::molecule::Molecule;
use std::collections::{HashSet, VecDeque};

/// A single ring in the molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    /// Atom indices forming the ring (in order around the ring).
    pub atoms: Vec<usize>,
    /// Bond indices forming the ring.
    pub bonds: Vec<usize>,
}

impl Ring {
    /// Returns the number of atoms in the ring.
    pub fn size(&self) -> usize {
        self.atoms.len()
    }

    /// Returns true if the given atom index is in this ring.
    pub fn contains_atom(&self, idx: usize) -> bool {
        self.atoms.contains(&idx)
    }

    /// Returns true if the given bond index is in this ring.
    pub fn contains_bond(&self, idx: usize) -> bool {
        self.bonds.contains(&idx)
    }
}

/// Compute the Smallest Set of Smallest Rings (SSSR).
///
/// Returns a minimum cycle basis (Horton's algorithm): candidate cycles are
/// built from shortest paths out of every ring atom, sorted by size, and
/// kept when linearly independent (over GF(2)) of the smaller ones. This is
/// guaranteed minimal, including for bridged systems such as cubane,
/// adamantane and norbornane. Ring order is deterministic.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::rings;
///
/// // Naphthalene: two fused 6-membered rings
/// let mut mol = Molecule::new("naphthalene");
/// for i in 0..10 {
///     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// // Ring 1: 0-1-2-3-4-5
/// for i in 0..5 {
///     mol.bonds.push(Bond::new(i, i + 1, BondOrder::Aromatic));
/// }
/// mol.bonds.push(Bond::new(5, 0, BondOrder::Aromatic));
/// // Ring 2: 3-4-6-7-8-9, shares edge 3-4
/// mol.bonds.push(Bond::new(4, 6, BondOrder::Aromatic));
/// mol.bonds.push(Bond::new(6, 7, BondOrder::Aromatic));
/// mol.bonds.push(Bond::new(7, 8, BondOrder::Aromatic));
/// mol.bonds.push(Bond::new(8, 9, BondOrder::Aromatic));
/// mol.bonds.push(Bond::new(9, 3, BondOrder::Aromatic));
///
/// let ring_set = rings::sssr(&mol);
/// assert_eq!(ring_set.len(), 2);
/// ```
pub fn sssr(mol: &Molecule) -> Vec<Ring> {
    let n = mol.atoms.len();
    let m = mol.bonds.len();
    if n == 0 || m == 0 {
        return vec![];
    }

    // Cyclomatic number: number of independent rings
    let num_components = connected_component_count(mol);
    if m + num_components <= n {
        return vec![];
    }
    let num_rings = m + num_components - n;

    let adj = AdjacencyList::from_molecule(mol);
    let ring_bond = crate::descriptors::topological::ring_bonds(mol);

    // Horton candidates: for each ring atom r and ring bond (x, y), the cycle
    // P(r, x) + (x, y) + P(y, r) using one BFS shortest-path tree per root.
    let mut candidates: Vec<Vec<usize>> = Vec::new();
    let mut seen: HashSet<Vec<usize>> = HashSet::new();
    let mut parent_bond = vec![usize::MAX; n];
    let mut dist = vec![usize::MAX; n];
    let mut on_path = vec![false; n];
    for root in 0..n {
        if !adj.neighbors(root).iter().any(|&(_, b)| ring_bond[b]) {
            continue;
        }
        // BFS over ring bonds only
        let mut reached = vec![root];
        dist[root] = 0;
        parent_bond[root] = usize::MAX;
        let mut queue = VecDeque::from([root]);
        while let Some(u) = queue.pop_front() {
            for &(v, b) in adj.neighbors(u) {
                if ring_bond[b] && dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    parent_bond[v] = b;
                    reached.push(v);
                    queue.push_back(v);
                }
            }
        }

        for (b, bond) in mol.bonds.iter().enumerate() {
            let (x, y) = (bond.atom1, bond.atom2);
            if !ring_bond[b]
                || x >= n
                || y >= n
                || dist[x] == usize::MAX
                || dist[y] == usize::MAX
                || parent_bond[x] == b
                || parent_bond[y] == b
            {
                continue;
            }
            let path_x = tree_path(mol, x, &parent_bond);
            let path_y = tree_path(mol, y, &parent_bond);
            // The two paths may only meet at the root
            for &bp in &path_x {
                let bnd = &mol.bonds[bp];
                on_path[bnd.atom1] = true;
                on_path[bnd.atom2] = true;
            }
            on_path[root] = false;
            let disjoint = !path_y.iter().any(|&bp| {
                let bnd = &mol.bonds[bp];
                (bnd.atom1 != root && on_path[bnd.atom1])
                    || (bnd.atom2 != root && on_path[bnd.atom2])
            }) && !(on_path[y] && y != root);
            for &bp in &path_x {
                let bnd = &mol.bonds[bp];
                on_path[bnd.atom1] = false;
                on_path[bnd.atom2] = false;
            }
            if !disjoint {
                continue;
            }
            let mut cycle: Vec<usize> = path_x.into_iter().chain(path_y).chain([b]).collect();
            cycle.sort_unstable();
            cycle.dedup();
            if cycle.len() >= 3 && seen.insert(cycle.clone()) {
                candidates.push(cycle);
            }
        }

        for v in reached {
            dist[v] = usize::MAX;
            parent_bond[v] = usize::MAX;
        }
    }

    candidates.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));

    // Greedy selection with incremental GF(2) elimination over bond bitsets
    let words = m.div_ceil(64);
    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new(); // (pivot bit, row)
    let mut rings = Vec::with_capacity(num_rings);
    for cycle in candidates {
        let mut row = vec![0u64; words];
        for &b in &cycle {
            row[b / 64] |= 1 << (b % 64);
        }
        for (pivot, basis_row) in &basis {
            if row[pivot / 64] >> (pivot % 64) & 1 == 1 {
                for (w, bw) in row.iter_mut().zip(basis_row) {
                    *w ^= bw;
                }
            }
        }
        let Some(pivot) = row
            .iter()
            .enumerate()
            .find(|(_, w)| **w != 0)
            .map(|(i, w)| i * 64 + w.trailing_zeros() as usize)
        else {
            continue; // dependent on smaller rings
        };
        basis.push((pivot, row));
        rings.push(ring_from_bonds(mol, &cycle));
        if rings.len() == num_rings {
            break;
        }
    }
    rings
}

/// Bonds on the BFS-tree path from `atom` back to the root.
fn tree_path(mol: &Molecule, mut atom: usize, parent_bond: &[usize]) -> Vec<usize> {
    let mut path = Vec::new();
    while parent_bond[atom] != usize::MAX {
        let b = parent_bond[atom];
        path.push(b);
        let bond = &mol.bonds[b];
        atom = if bond.atom1 == atom {
            bond.atom2
        } else {
            bond.atom1
        };
    }
    path
}

/// Builds a [`Ring`] with atoms and bonds in walking order from a cycle's bonds.
fn ring_from_bonds(mol: &Molecule, cycle: &[usize]) -> Ring {
    let start_bond = cycle[0];
    let start = mol.bonds[start_bond].atom1;
    let mut atoms = vec![start];
    let mut bonds = vec![start_bond];
    let mut current = mol.bonds[start_bond].atom2;
    while current != start && bonds.len() < cycle.len() {
        atoms.push(current);
        let next = cycle
            .iter()
            .copied()
            .find(|&b| !bonds.contains(&b) && mol.bonds[b].contains_atom(current));
        let Some(b) = next else { break };
        bonds.push(b);
        current = mol.bonds[b].other_atom(current).unwrap_or(start);
    }
    Ring { atoms, bonds }
}

/// Count connected components using union-find.
fn connected_component_count(mol: &Molecule) -> usize {
    if mol.atoms.is_empty() {
        return 0;
    }
    let n = mol.atoms.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    for bond in &mol.bonds {
        if bond.atom1 < n && bond.atom2 < n {
            let pa = find(&mut parent, bond.atom1);
            let pb = find(&mut parent, bond.atom2);
            if pa != pb {
                parent[pa] = pb;
            }
        }
    }

    let mut roots = HashSet::new();
    for i in 0..n {
        roots.insert(find(&mut parent, i));
    }
    roots.len()
}

/// Get the ring sizes for a specific atom.
///
/// Returns a vector of ring sizes that the atom participates in.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::rings;
///
/// let mut mol = Molecule::new("benzene");
/// for i in 0..6 {
///     mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     mol.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
///
/// let sizes = rings::ring_sizes(&mol, 0);
/// assert_eq!(sizes, vec![6]);
/// ```
pub fn ring_sizes(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let rings = sssr(mol);
    rings
        .iter()
        .filter(|r| r.contains_atom(atom_idx))
        .map(|r| r.size())
        .collect()
}

/// Get the smallest ring size for a specific atom.
///
/// Returns `None` if the atom is not in any ring.
pub fn smallest_ring_size(mol: &Molecule, atom_idx: usize) -> Option<usize> {
    ring_sizes(mol, atom_idx).into_iter().min()
}

/// Check if an atom is in a ring of a specific size.
pub fn is_in_ring_of_size(mol: &Molecule, atom_idx: usize, size: usize) -> bool {
    ring_sizes(mol, atom_idx).contains(&size)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    fn make_benzene() -> Molecule {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            let angle = std::f64::consts::PI * 2.0 * i as f64 / 6.0;
            mol.atoms
                .push(Atom::new(i, "C", angle.cos(), angle.sin(), 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }
        mol
    }

    fn make_cyclopentane() -> Molecule {
        let mut mol = Molecule::new("cyclopentane");
        for i in 0..5 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..5 {
            mol.bonds.push(Bond::new(i, (i + 1) % 5, BondOrder::Single));
        }
        mol
    }

    fn make_naphthalene() -> Molecule {
        let mut mol = Molecule::new("naphthalene");
        for i in 0..10 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        // Ring 1: 0-1-2-3-4-5
        mol.bonds.push(Bond::new(0, 1, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(3, 4, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(4, 5, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(5, 0, BondOrder::Aromatic));
        // Ring 2: 4-6-7-8-9-3 (shares edge 3-4)
        mol.bonds.push(Bond::new(4, 6, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(6, 7, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(7, 8, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(8, 9, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(9, 3, BondOrder::Aromatic));
        mol
    }

    #[test]
    fn test_sssr_benzene() {
        let mol = make_benzene();
        let rings = sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].size(), 6);
    }

    #[test]
    fn test_sssr_cyclopentane() {
        let mol = make_cyclopentane();
        let rings = sssr(&mol);
        assert_eq!(rings.len(), 1);
        assert_eq!(rings[0].size(), 5);
    }

    #[test]
    fn test_sssr_naphthalene() {
        let mol = make_naphthalene();
        let rings = sssr(&mol);
        assert_eq!(rings.len(), 2);
        // Both should be 6-membered
        for ring in &rings {
            assert_eq!(ring.size(), 6);
        }
    }

    #[test]
    fn test_sssr_propane() {
        let mut mol = Molecule::new("propane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 3.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));

        let rings = sssr(&mol);
        assert_eq!(rings.len(), 0);
    }

    #[test]
    fn test_sssr_empty() {
        let mol = Molecule::new("empty");
        let rings = sssr(&mol);
        assert_eq!(rings.len(), 0);
    }

    #[test]
    fn test_ring_sizes_benzene() {
        let mol = make_benzene();
        for i in 0..6 {
            let sizes = ring_sizes(&mol, i);
            assert_eq!(sizes, vec![6]);
        }
    }

    #[test]
    fn test_smallest_ring_size() {
        let mol = make_benzene();
        assert_eq!(smallest_ring_size(&mol, 0), Some(6));
    }

    #[test]
    fn test_smallest_ring_size_acyclic() {
        let mut mol = Molecule::new("ethane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        assert_eq!(smallest_ring_size(&mol, 0), None);
    }

    #[test]
    fn test_is_in_ring_of_size() {
        let mol = make_benzene();
        assert!(is_in_ring_of_size(&mol, 0, 6));
        assert!(!is_in_ring_of_size(&mol, 0, 5));
    }

    #[test]
    fn test_ring_contains_atom() {
        let ring = Ring {
            atoms: vec![0, 1, 2, 3, 4, 5],
            bonds: vec![0, 1, 2, 3, 4, 5],
        };
        assert!(ring.contains_atom(3));
        assert!(!ring.contains_atom(7));
    }

    #[test]
    fn test_ring_contains_bond() {
        let ring = Ring {
            atoms: vec![0, 1, 2],
            bonds: vec![0, 1, 2],
        };
        assert!(ring.contains_bond(1));
        assert!(!ring.contains_bond(5));
    }

    #[test]
    fn test_cubane_rings() {
        // Cubane: 8 atoms, 12 bonds, 5 independent rings
        let mut mol = Molecule::new("cubane");
        for i in 0..8 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        // Bottom face: 0-1-2-3
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(3, 0, BondOrder::Single));
        // Top face: 4-5-6-7
        mol.bonds.push(Bond::new(4, 5, BondOrder::Single));
        mol.bonds.push(Bond::new(5, 6, BondOrder::Single));
        mol.bonds.push(Bond::new(6, 7, BondOrder::Single));
        mol.bonds.push(Bond::new(7, 4, BondOrder::Single));
        // Vertical: 0-4, 1-5, 2-6, 3-7
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 5, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 6, BondOrder::Single));
        mol.bonds.push(Bond::new(3, 7, BondOrder::Single));

        let rings = sssr(&mol);
        // Cubane: 12 - 8 + 1 = 5 independent rings
        assert_eq!(rings.len(), 5);
    }
}
