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
/// Uses a shortest-path based approach: for each edge not in the spanning tree,
/// find the shortest cycle through that edge using BFS.
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

    let adj = AdjacencyList::from_molecule(mol);
    let num_components = connected_component_count(mol);

    // Cyclomatic number: number of independent rings
    let num_rings = if m + num_components > n {
        m - n + num_components
    } else {
        return vec![];
    };

    if num_rings == 0 {
        return vec![];
    }

    // Build spanning tree using BFS, collect non-tree edges
    let mut in_tree = vec![false; m];
    let mut visited = vec![false; n];
    let mut non_tree_edges: Vec<usize> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut queue = VecDeque::new();
        visited[start] = true;
        queue.push_back(start);

        while let Some(u) = queue.pop_front() {
            for &(v, bond_idx) in adj.neighbors(u) {
                if !visited[v] {
                    visited[v] = true;
                    in_tree[bond_idx] = true;
                    queue.push_back(v);
                }
            }
        }
    }

    // Non-tree edges form cycles
    for (bond_idx, &is_tree) in in_tree.iter().enumerate() {
        if !is_tree {
            non_tree_edges.push(bond_idx);
        }
    }

    // For each non-tree edge, find the shortest cycle containing it
    let mut candidate_rings: Vec<(usize, Vec<usize>)> = Vec::new(); // (size, atom_path)

    for &edge_idx in &non_tree_edges {
        let bond = &mol.bonds[edge_idx];
        // Find shortest path from bond.atom1 to bond.atom2 NOT using this edge
        if let Some(path) = shortest_path_avoiding_edge(&adj, bond.atom1, bond.atom2, edge_idx, n) {
            candidate_rings.push((path.len(), path));
        }
    }

    // Sort by ring size (smallest first) for SSSR selection
    candidate_rings.sort_by_key(|(size, _)| *size);

    // Select linearly independent rings using bond-set independence
    let mut selected_rings: Vec<Vec<usize>> = Vec::new();
    let mut selected_bond_sets: Vec<HashSet<usize>> = Vec::new();

    for (_size, atom_path) in &candidate_rings {
        if selected_rings.len() >= num_rings {
            break;
        }

        let bond_set = path_to_bond_set(atom_path, mol, &adj);

        if is_linearly_independent(&bond_set, &selected_bond_sets) {
            selected_rings.push(atom_path.clone());
            selected_bond_sets.push(bond_set);
        }
    }

    // If we don't have enough from non-tree edges, also try BFS-found rings
    if selected_rings.len() < num_rings {
        // Find additional rings by shortest path between endpoints of remaining non-tree edges
        // using different spanning trees or by enumerating small cycles
        let extra = find_additional_rings(
            mol,
            &adj,
            num_rings - selected_rings.len(),
            &selected_bond_sets,
        );
        for (path, bond_set) in extra {
            if selected_rings.len() >= num_rings {
                break;
            }
            if is_linearly_independent(&bond_set, &selected_bond_sets) {
                selected_rings.push(path);
                selected_bond_sets.push(bond_set);
            }
        }
    }

    // Convert to Ring structs
    selected_rings
        .iter()
        .zip(selected_bond_sets.iter())
        .map(|(atoms, bonds)| Ring {
            atoms: atoms.clone(),
            bonds: bonds.iter().copied().collect(),
        })
        .collect()
}

/// Find additional rings using BFS from every atom to find small cycles.
fn find_additional_rings(
    mol: &Molecule,
    adj: &AdjacencyList,
    needed: usize,
    existing_bond_sets: &[HashSet<usize>],
) -> Vec<(Vec<usize>, HashSet<usize>)> {
    let n = adj.num_atoms();
    let mut result = Vec::new();
    let mut found_bond_sets: Vec<HashSet<usize>> = existing_bond_sets.to_vec();

    // For each atom, do BFS to find short cycles
    for start in 0..n {
        if result.len() >= needed {
            break;
        }

        let cycles = find_short_cycles_from(start, adj, n);
        for cycle in cycles {
            if result.len() >= needed {
                break;
            }
            let bond_set = path_to_bond_set(&cycle, mol, adj);
            if is_linearly_independent(&bond_set, &found_bond_sets) {
                found_bond_sets.push(bond_set.clone());
                result.push((cycle, bond_set));
            }
        }
    }

    result
}

/// Find short cycles starting from a given atom using BFS.
fn find_short_cycles_from(start: usize, adj: &AdjacencyList, n: usize) -> Vec<Vec<usize>> {
    let mut cycles = Vec::new();
    let mut dist = vec![usize::MAX; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    dist[start] = 0;
    queue.push_back(start);

    while let Some(u) = queue.pop_front() {
        for &(v, _bond_idx) in adj.neighbors(u) {
            if dist[v] == usize::MAX {
                dist[v] = dist[u] + 1;
                parent[v] = u;
                queue.push_back(v);
            } else if v != parent[u] && dist[v] >= dist[u] {
                // Found a cycle: trace paths from u and v back to their LCA
                let path_u = trace_to_root(u, &parent);
                let path_v = trace_to_root(v, &parent);

                // Find where the paths diverge (LCA)
                let cycle = merge_paths_to_cycle(&path_u, &path_v);
                if cycle.len() >= 3 {
                    // Deduplicate
                    let normalized = normalize_cycle(&cycle);
                    cycles.push(normalized);
                }
            }
        }
    }

    cycles
}

/// Trace path from node to root (start of BFS).
fn trace_to_root(node: usize, parent: &[usize]) -> Vec<usize> {
    let mut path = vec![node];
    let mut current = node;
    while parent[current] != usize::MAX {
        current = parent[current];
        path.push(current);
    }
    path
}

/// Merge two paths from their common root to form a cycle.
fn merge_paths_to_cycle(path_u: &[usize], path_v: &[usize]) -> Vec<usize> {
    // Both paths end at the same root (BFS start)
    // Find the longest common suffix (from root)
    let mut i = path_u.len();
    let mut j = path_v.len();

    // Walk from root to find divergence point
    while i > 0 && j > 0 && path_u[i - 1] == path_v[j - 1] {
        i -= 1;
        j -= 1;
    }

    // The divergence point
    let lca_idx_u = i; // path_u[i..] is common, path_u[..=i] is the unique part (inclusive of LCA)
    let lca_idx_v = j;

    // Cycle = path_u[0..=lca_idx_u] + path_v[0..lca_idx_v] reversed
    let mut cycle: Vec<usize> = path_u[..=lca_idx_u].to_vec();
    for k in (0..lca_idx_v).rev() {
        cycle.push(path_v[k]);
    }

    cycle
}

/// Normalize a cycle for deduplication (rotate to start at smallest index).
fn normalize_cycle(cycle: &[usize]) -> Vec<usize> {
    if cycle.is_empty() {
        return vec![];
    }
    let min_pos = cycle.iter().enumerate().min_by_key(|&(_, v)| *v).unwrap().0;
    let mut normalized = Vec::with_capacity(cycle.len());
    for i in 0..cycle.len() {
        normalized.push(cycle[(min_pos + i) % cycle.len()]);
    }
    normalized
}

/// Find shortest path from `src` to `dst` in the graph, avoiding `excluded_edge`.
fn shortest_path_avoiding_edge(
    adj: &AdjacencyList,
    src: usize,
    dst: usize,
    excluded_edge: usize,
    n: usize,
) -> Option<Vec<usize>> {
    let mut dist = vec![usize::MAX; n];
    let mut parent = vec![usize::MAX; n];
    let mut queue = VecDeque::new();

    dist[src] = 0;
    queue.push_back(src);

    while let Some(u) = queue.pop_front() {
        if u == dst {
            // Reconstruct path
            let mut path = vec![dst];
            let mut current = dst;
            while current != src {
                current = parent[current];
                path.push(current);
            }
            path.reverse();
            return Some(path);
        }

        for &(v, bond_idx) in adj.neighbors(u) {
            if bond_idx == excluded_edge {
                continue;
            }
            if dist[v] == usize::MAX {
                dist[v] = dist[u] + 1;
                parent[v] = u;
                queue.push_back(v);
            }
        }
    }

    None
}

/// Convert an atom path (cycle) to a set of bond indices.
fn path_to_bond_set(atom_path: &[usize], _mol: &Molecule, adj: &AdjacencyList) -> HashSet<usize> {
    let mut bonds = HashSet::new();
    let len = atom_path.len();
    for i in 0..len {
        let a = atom_path[i];
        let b = atom_path[(i + 1) % len];
        // Find bond between a and b
        for &(neighbor, bond_idx) in adj.neighbors(a) {
            if neighbor == b {
                bonds.insert(bond_idx);
                break;
            }
        }
    }
    bonds
}

/// Check if a bond set is linearly independent of existing bond sets.
///
/// Uses GF(2) vector space: a ring is independent if its bond set cannot
/// be expressed as the symmetric difference of a subset of existing rings.
fn is_linearly_independent(new_bonds: &HashSet<usize>, existing: &[HashSet<usize>]) -> bool {
    if existing.is_empty() {
        return true;
    }

    // Use Gaussian elimination over GF(2)
    // Convert bond sets to bit vectors
    let max_bond = new_bonds
        .iter()
        .copied()
        .chain(existing.iter().flat_map(|s| s.iter().copied()))
        .max()
        .unwrap_or(0)
        + 1;

    let to_bitvec = |s: &HashSet<usize>| -> Vec<bool> {
        let mut bv = vec![false; max_bond];
        for &b in s {
            bv[b] = true;
        }
        bv
    };

    // Build matrix from existing + new
    let mut matrix: Vec<Vec<bool>> = existing.iter().map(&to_bitvec).collect();
    matrix.push(to_bitvec(new_bonds));

    // Gaussian elimination over GF(2)
    let rows = matrix.len();
    let cols = max_bond;
    let mut pivot_row = 0;

    for col in 0..cols {
        // Find pivot
        let mut found = false;
        for row in pivot_row..rows {
            if matrix[row][col] {
                matrix.swap(pivot_row, row);
                found = true;
                break;
            }
        }
        if !found {
            continue;
        }

        // Eliminate
        let pivot_vals: Vec<bool> = matrix[pivot_row].clone();
        for (row, row_data) in matrix.iter_mut().enumerate() {
            if row != pivot_row && row_data[col] {
                for (c, &pv) in pivot_vals.iter().enumerate() {
                    row_data[c] ^= pv;
                }
            }
        }
        pivot_row += 1;
    }

    // The new ring (last row) is independent if it wasn't eliminated to zero
    matrix[rows - 1].iter().any(|&b| b)
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
