//! Topological molecular descriptors.
//!
//! This module provides functions for calculating topological properties
//! such as ring count, ring membership, and rotatable bond count.

use crate::bond::BondOrder;
use crate::molecule::Molecule;

/// Count the number of rings in the molecule.
///
/// Uses the Euler characteristic formula: rings = bonds - atoms + components.
/// This gives the number of independent cycles (cyclomatic number).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder, descriptors};
///
/// // Benzene has 1 ring
/// let mut benzene = Molecule::new("benzene");
/// for i in 0..6 {
///     benzene.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     benzene.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
/// assert_eq!(descriptors::ring_count(&benzene), 1);
/// ```
pub fn ring_count(mol: &Molecule) -> usize {
    if mol.atoms.is_empty() {
        return 0;
    }

    let num_atoms = mol.atoms.len();
    let num_bonds = mol.bonds.len();
    let num_components = connected_component_count(mol);

    // Euler formula for cyclomatic number: E - V + C
    // where E = edges (bonds), V = vertices (atoms), C = connected components
    if num_bonds + num_components > num_atoms {
        num_bonds - num_atoms + num_components
    } else {
        0
    }
}

/// Get ring membership for all atoms.
///
/// Returns a vector of booleans where `result[i]` is `true` if atom `i` is in a ring.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder, descriptors};
///
/// // Benzene: all 6 carbons are in the ring
/// let mut benzene = Molecule::new("benzene");
/// for i in 0..6 {
///     benzene.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     benzene.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
///
/// let ring_atoms = descriptors::ring_atoms(&benzene);
/// assert!(ring_atoms.iter().all(|&in_ring| in_ring));
/// ```
pub fn ring_atoms(mol: &Molecule) -> Vec<bool> {
    let (atom_in_ring, _) = compute_ring_membership(mol);
    atom_in_ring
}

/// Get ring membership for all bonds.
///
/// Returns a vector of booleans where `result[i]` is `true` if bond `i` is in a ring.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder, descriptors};
///
/// // Benzene: all 6 bonds are in the ring
/// let mut benzene = Molecule::new("benzene");
/// for i in 0..6 {
///     benzene.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     benzene.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
///
/// let ring_bonds = descriptors::ring_bonds(&benzene);
/// assert!(ring_bonds.iter().all(|&in_ring| in_ring));
/// ```
pub fn ring_bonds(mol: &Molecule) -> Vec<bool> {
    let (_, bond_in_ring) = compute_ring_membership(mol);
    bond_in_ring
}

/// Count rotatable bonds.
///
/// A bond is considered rotatable if:
/// - It is a single bond (BondOrder::Single)
/// - It is not in a ring
/// - Neither atom is hydrogen
/// - Both atoms are non-terminal (each is bonded to at least one other heavy atom)
///
/// This follows the common definition: "a rotatable bond is any single non-ring bond,
/// bounded to nonterminal heavy atoms."
///
/// Terminal atoms are those bonded only to hydrogens (and this one heavy atom bond).
/// For example, methyl groups (-CH3), amino groups (-NH2), and hydroxyl groups (-OH)
/// are terminal.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder, descriptors};
///
/// // Benzene has 0 rotatable bonds (all in ring, and aromatic not single)
/// let mut benzene = Molecule::new("benzene");
/// for i in 0..6 {
///     benzene.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
/// }
/// for i in 0..6 {
///     benzene.bonds.push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
/// }
/// assert_eq!(descriptors::rotatable_bond_count(&benzene), 0);
/// ```
pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    if mol.bonds.is_empty() {
        return 0;
    }

    let (_, bond_in_ring) = compute_ring_membership(mol);
    let heavy_degree = compute_heavy_atom_degrees(mol);

    mol.bonds
        .iter()
        .enumerate()
        .filter(|(i, bond)| {
            // Must be a single bond
            if bond.order != BondOrder::Single {
                return false;
            }

            // Must not be in a ring
            if bond_in_ring[*i] {
                return false;
            }

            // Neither atom can be hydrogen
            let atom1_elem = mol.atoms.get(bond.atom1).map(|a| a.element.as_str());
            let atom2_elem = mol.atoms.get(bond.atom2).map(|a| a.element.as_str());
            if is_hydrogen(atom1_elem) || is_hydrogen(atom2_elem) {
                return false;
            }

            // Both atoms must be non-terminal (heavy_degree > 1)
            // heavy_degree counts bonds to other heavy atoms
            let hdeg1 = heavy_degree.get(bond.atom1).copied().unwrap_or(0);
            let hdeg2 = heavy_degree.get(bond.atom2).copied().unwrap_or(0);
            if hdeg1 <= 1 || hdeg2 <= 1 {
                return false;
            }

            true
        })
        .count()
}

/// Count connected components using union-find.
fn connected_component_count(mol: &Molecule) -> usize {
    if mol.atoms.is_empty() {
        return 0;
    }

    let n = mol.atoms.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0usize; n];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], rank: &mut [usize], x: usize, y: usize) {
        let px = find(parent, x);
        let py = find(parent, y);
        if px != py {
            if rank[px] < rank[py] {
                parent[px] = py;
            } else if rank[px] > rank[py] {
                parent[py] = px;
            } else {
                parent[py] = px;
                rank[px] += 1;
            }
        }
    }

    for bond in &mol.bonds {
        if bond.atom1 < n && bond.atom2 < n {
            union(&mut parent, &mut rank, bond.atom1, bond.atom2);
        }
    }

    // Count unique roots
    let mut roots = std::collections::HashSet::new();
    for i in 0..n {
        roots.insert(find(&mut parent, i));
    }
    roots.len()
}

/// Compute ring membership for atoms and bonds.
///
/// Uses DFS to find all cycles and marks atoms/bonds that are part of any cycle.
fn compute_ring_membership(mol: &Molecule) -> (Vec<bool>, Vec<bool>) {
    let n = mol.atoms.len();
    let m = mol.bonds.len();

    if n == 0 {
        return (vec![], vec![]);
    }

    let mut atom_in_ring = vec![false; n];
    let mut bond_in_ring = vec![false; m];

    // Build adjacency list with bond indices
    let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n]; // (neighbor, bond_index)
    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        if bond.atom1 < n && bond.atom2 < n {
            adj[bond.atom1].push((bond.atom2, bond_idx));
            adj[bond.atom2].push((bond.atom1, bond_idx));
        }
    }

    // For each connected component, do DFS to find back edges
    let mut visited = vec![false; n];
    let mut in_stack = vec![false; n];
    let mut parent = vec![usize::MAX; n];
    let mut parent_bond = vec![usize::MAX; n];

    for start in 0..n {
        if visited[start] {
            continue;
        }

        // DFS with explicit stack: (node, neighbor_index, entering)
        // entering=true means we're visiting this node for the first time
        let mut stack: Vec<(usize, usize, bool)> = vec![(start, 0, true)];

        while let Some((node, next_neighbor_idx, entering)) = stack.pop() {
            if entering {
                if visited[node] {
                    continue;
                }
                visited[node] = true;
                in_stack[node] = true;
                // Re-push with entering=false to handle backtracking
                stack.push((node, 0, false));
            } else {
                // Process neighbors
                if next_neighbor_idx < adj[node].len() {
                    let (neighbor, bond_idx) = adj[node][next_neighbor_idx];
                    // Push current node back with next neighbor index
                    stack.push((node, next_neighbor_idx + 1, false));

                    // Skip the edge we came from
                    if parent_bond[node] == bond_idx {
                        continue;
                    }

                    if !visited[neighbor] {
                        // Tree edge - go deeper
                        parent[neighbor] = node;
                        parent_bond[neighbor] = bond_idx;
                        stack.push((neighbor, 0, true));
                    } else if in_stack[neighbor] {
                        // Back edge found - mark the cycle
                        mark_cycle_path(
                            node,
                            neighbor,
                            bond_idx,
                            &parent,
                            &parent_bond,
                            &mut atom_in_ring,
                            &mut bond_in_ring,
                        );
                    }
                } else {
                    // Done with this node, remove from stack
                    in_stack[node] = false;
                }
            }
        }
    }

    (atom_in_ring, bond_in_ring)
}

/// Mark atoms and bonds in a cycle found via back edge from u to v.
fn mark_cycle_path(
    u: usize,
    v: usize,
    back_edge_bond: usize,
    parent: &[usize],
    parent_bond: &[usize],
    atom_in_ring: &mut [bool],
    bond_in_ring: &mut [bool],
) {
    // Mark the back edge itself
    bond_in_ring[back_edge_bond] = true;

    // Trace from u back to v through parent pointers
    let mut current = u;
    while current != v && current != usize::MAX {
        atom_in_ring[current] = true;
        if parent_bond[current] != usize::MAX {
            bond_in_ring[parent_bond[current]] = true;
        }
        current = parent[current];
    }
    // Mark v
    if current == v {
        atom_in_ring[v] = true;
    }
}

/// Compute degree (number of bonds) for each atom.
#[allow(dead_code)]
fn compute_degrees(mol: &Molecule) -> Vec<usize> {
    let n = mol.atoms.len();
    let mut degrees = vec![0usize; n];
    for bond in &mol.bonds {
        if bond.atom1 < n {
            degrees[bond.atom1] += 1;
        }
        if bond.atom2 < n {
            degrees[bond.atom2] += 1;
        }
    }
    degrees
}

/// Compute heavy atom degree (number of bonds to non-hydrogen atoms) for each atom.
fn compute_heavy_atom_degrees(mol: &Molecule) -> Vec<usize> {
    let n = mol.atoms.len();
    let mut degrees = vec![0usize; n];
    for bond in &mol.bonds {
        if bond.atom1 < n && bond.atom2 < n {
            let elem1 = mol.atoms.get(bond.atom1).map(|a| a.element.as_str());
            let elem2 = mol.atoms.get(bond.atom2).map(|a| a.element.as_str());

            // If atom2 is not hydrogen, increment atom1's heavy degree
            if !is_hydrogen(elem2) {
                degrees[bond.atom1] += 1;
            }
            // If atom1 is not hydrogen, increment atom2's heavy degree
            if !is_hydrogen(elem1) {
                degrees[bond.atom2] += 1;
            }
        }
    }
    degrees
}

/// Check if an element is hydrogen.
fn is_hydrogen(element: Option<&str>) -> bool {
    match element {
        Some(e) => {
            let e = e.trim().to_uppercase();
            e == "H" || e == "D" || e == "T"
        }
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::Bond;

    fn make_benzene() -> Molecule {
        let mut mol = Molecule::new("benzene");
        // Add 6 carbon atoms
        for i in 0..6 {
            let angle = std::f64::consts::PI * 2.0 * i as f64 / 6.0;
            mol.atoms
                .push(Atom::new(i, "C", angle.cos(), angle.sin(), 0.0));
        }
        // Add 6 aromatic bonds in a ring
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }
        mol
    }

    fn make_propane() -> Molecule {
        let mut mol = Molecule::new("propane");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 3.0, 0.0, 0.0));
        // Add hydrogens
        for i in 3..11 {
            mol.atoms.push(Atom::new(i, "H", 0.0, 0.0, 0.0));
        }
        // C-C bonds
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        // C-H bonds (simplified)
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 5, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 6, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 7, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 8, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 9, BondOrder::Single));
        mol.bonds.push(Bond::new(2, 10, BondOrder::Single));
        mol
    }

    fn make_naphthalene() -> Molecule {
        // Two fused benzene rings
        let mut mol = Molecule::new("naphthalene");
        for i in 0..10 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        // First ring: 0-1-2-3-4-5
        mol.bonds.push(Bond::new(0, 1, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(2, 3, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(3, 4, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(4, 5, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(5, 0, BondOrder::Aromatic));
        // Second ring: 3-4-6-7-8-9 (shares 3-4)
        mol.bonds.push(Bond::new(4, 6, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(6, 7, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(7, 8, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(8, 9, BondOrder::Aromatic));
        mol.bonds.push(Bond::new(9, 3, BondOrder::Aromatic));
        mol
    }

    #[test]
    fn test_ring_count_benzene() {
        let mol = make_benzene();
        assert_eq!(ring_count(&mol), 1);
    }

    #[test]
    fn test_ring_count_propane() {
        let mol = make_propane();
        assert_eq!(ring_count(&mol), 0);
    }

    #[test]
    fn test_ring_count_naphthalene() {
        let mol = make_naphthalene();
        assert_eq!(ring_count(&mol), 2);
    }

    #[test]
    fn test_ring_count_empty() {
        let mol = Molecule::new("empty");
        assert_eq!(ring_count(&mol), 0);
    }

    #[test]
    fn test_ring_count_single_atom() {
        let mut mol = Molecule::new("single");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        assert_eq!(ring_count(&mol), 0);
    }

    #[test]
    fn test_ring_atoms_benzene() {
        let mol = make_benzene();
        let in_ring = ring_atoms(&mol);
        assert_eq!(in_ring.len(), 6);
        assert!(in_ring.iter().all(|&r| r), "All atoms in benzene should be in ring");
    }

    #[test]
    fn test_ring_atoms_propane() {
        let mol = make_propane();
        let in_ring = ring_atoms(&mol);
        assert!(
            in_ring.iter().all(|&r| !r),
            "No atoms in propane should be in ring"
        );
    }

    #[test]
    fn test_ring_bonds_benzene() {
        let mol = make_benzene();
        let in_ring = ring_bonds(&mol);
        assert_eq!(in_ring.len(), 6);
        assert!(in_ring.iter().all(|&r| r), "All bonds in benzene should be in ring");
    }

    #[test]
    fn test_ring_bonds_propane() {
        let mol = make_propane();
        let in_ring = ring_bonds(&mol);
        assert!(
            in_ring.iter().all(|&r| !r),
            "No bonds in propane should be in ring"
        );
    }

    #[test]
    fn test_rotatable_bond_count_benzene() {
        let mol = make_benzene();
        // All aromatic bonds, not single, so 0 rotatable
        assert_eq!(rotatable_bond_count(&mol), 0);
    }

    #[test]
    fn test_rotatable_bond_count_propane() {
        let mol = make_propane();
        // By RDKit definition, rotatable bonds require both atoms to be non-terminal
        // (heavy_degree > 1). In propane:
        // - C1 (methyl) has heavy_degree = 1 (only bonded to C2)
        // - C2 (middle) has heavy_degree = 2 (bonded to C1 and C3)
        // - C3 (methyl) has heavy_degree = 1 (only bonded to C2)
        // Both C1-C2 and C2-C3 bonds have one terminal atom, so neither is rotatable
        assert_eq!(rotatable_bond_count(&mol), 0);
    }

    #[test]
    fn test_rotatable_bond_count_empty() {
        let mol = Molecule::new("empty");
        assert_eq!(rotatable_bond_count(&mol), 0);
    }

    #[test]
    fn test_rotatable_bond_count_ethane_skeleton() {
        // Ethane without hydrogens - just C-C
        let mut mol = Molecule::new("ethane_skeleton");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        // Both atoms are terminal (degree 1), so not rotatable
        assert_eq!(rotatable_bond_count(&mol), 0);
    }

    #[test]
    fn test_connected_components_single() {
        let mol = make_benzene();
        assert_eq!(connected_component_count(&mol), 1);
    }

    #[test]
    fn test_connected_components_multiple() {
        let mut mol = Molecule::new("two_fragments");
        // Fragment 1: C-C
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        // Fragment 2: O (isolated)
        mol.atoms.push(Atom::new(2, "O", 5.0, 0.0, 0.0));
        assert_eq!(connected_component_count(&mol), 2);
    }

    #[test]
    fn test_connected_components_empty() {
        let mol = Molecule::new("empty");
        assert_eq!(connected_component_count(&mol), 0);
    }
}
