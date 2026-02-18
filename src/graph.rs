//! Graph adjacency infrastructure for molecular structures.
//!
//! Provides pre-computed adjacency lists and degree vectors for efficient
//! graph-based algorithms (ring perception, featurization, fingerprints).
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::graph::AdjacencyList;
//!
//! let mut mol = Molecule::new("methane");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! let adj = AdjacencyList::from_molecule(&mol);
//! assert_eq!(adj.degree(0), 2);
//! assert_eq!(adj.neighbors(0), &[(1, 0), (2, 1)]);
//! ```

use crate::molecule::Molecule;

/// Pre-computed adjacency list for a molecule's bond graph.
///
/// Stores neighbor lists and degree vectors for O(1) access.
/// Each neighbor entry is `(atom_index, bond_index)`.
#[derive(Debug, Clone)]
pub struct AdjacencyList {
    /// For each atom, a list of (neighbor_atom_index, bond_index) pairs.
    adj: Vec<Vec<(usize, usize)>>,
    /// Degree (total bond count) for each atom.
    degrees: Vec<usize>,
    /// Heavy atom degree (bonds to non-H atoms) for each atom.
    heavy_degrees: Vec<usize>,
    /// Number of atoms.
    num_atoms: usize,
}

impl AdjacencyList {
    /// Build an adjacency list from a molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use sdfrust::{Molecule, Atom, Bond, BondOrder};
    /// use sdfrust::graph::AdjacencyList;
    ///
    /// let mut mol = Molecule::new("water");
    /// mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
    /// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
    /// mol.atoms.push(Atom::new(2, "H", -0.3, 0.9, 0.0));
    /// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
    /// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
    ///
    /// let adj = AdjacencyList::from_molecule(&mol);
    /// assert_eq!(adj.num_atoms(), 3);
    /// assert_eq!(adj.degree(0), 2);
    /// assert_eq!(adj.degree(1), 1);
    /// ```
    pub fn from_molecule(mol: &Molecule) -> Self {
        let n = mol.atoms.len();
        let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n];
        let mut degrees = vec![0usize; n];
        let mut heavy_degrees = vec![0usize; n];

        for (bond_idx, bond) in mol.bonds.iter().enumerate() {
            if bond.atom1 < n && bond.atom2 < n {
                adj[bond.atom1].push((bond.atom2, bond_idx));
                adj[bond.atom2].push((bond.atom1, bond_idx));
                degrees[bond.atom1] += 1;
                degrees[bond.atom2] += 1;

                let elem1 = mol.atoms[bond.atom1].element.as_str();
                let elem2 = mol.atoms[bond.atom2].element.as_str();
                if !is_hydrogen(elem2) {
                    heavy_degrees[bond.atom1] += 1;
                }
                if !is_hydrogen(elem1) {
                    heavy_degrees[bond.atom2] += 1;
                }
            }
        }

        Self {
            adj,
            degrees,
            heavy_degrees,
            num_atoms: n,
        }
    }

    /// Returns the number of atoms in the graph.
    pub fn num_atoms(&self) -> usize {
        self.num_atoms
    }

    /// Returns the neighbors of atom `idx` as `(neighbor_atom, bond_index)` pairs.
    ///
    /// Returns an empty slice if the index is out of bounds.
    pub fn neighbors(&self, idx: usize) -> &[(usize, usize)] {
        self.adj.get(idx).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Returns just the neighbor atom indices for atom `idx`.
    pub fn neighbor_atoms(&self, idx: usize) -> Vec<usize> {
        self.neighbors(idx).iter().map(|(a, _)| *a).collect()
    }

    /// Returns the degree (number of bonds) for atom `idx`.
    ///
    /// Returns 0 if the index is out of bounds.
    pub fn degree(&self, idx: usize) -> usize {
        self.degrees.get(idx).copied().unwrap_or(0)
    }

    /// Returns the heavy atom degree (bonds to non-H atoms) for atom `idx`.
    ///
    /// Returns 0 if the index is out of bounds.
    pub fn heavy_degree(&self, idx: usize) -> usize {
        self.heavy_degrees.get(idx).copied().unwrap_or(0)
    }

    /// Returns a slice of all degrees.
    pub fn degrees(&self) -> &[usize] {
        &self.degrees
    }

    /// Returns a slice of all heavy atom degrees.
    pub fn heavy_degrees(&self) -> &[usize] {
        &self.heavy_degrees
    }

    /// Returns the bond indices for bonds involving atom `idx`.
    pub fn bond_indices(&self, idx: usize) -> Vec<usize> {
        self.neighbors(idx).iter().map(|(_, b)| *b).collect()
    }
}

/// Check if an element symbol is hydrogen (H, D, or T).
pub(crate) fn is_hydrogen(element: &str) -> bool {
    let e = element.trim();
    e.eq_ignore_ascii_case("H") || e.eq_ignore_ascii_case("D") || e.eq_ignore_ascii_case("T")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    #[test]
    fn test_adjacency_list_methane() {
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

        let adj = AdjacencyList::from_molecule(&mol);
        assert_eq!(adj.num_atoms(), 5);
        assert_eq!(adj.degree(0), 4);
        assert_eq!(adj.heavy_degree(0), 0); // All neighbors are H
        assert_eq!(adj.degree(1), 1);
        assert_eq!(adj.heavy_degree(1), 1); // C is a heavy atom
    }

    #[test]
    fn test_adjacency_list_benzene() {
        let mut mol = Molecule::new("benzene");
        for i in 0..6 {
            mol.atoms.push(Atom::new(i, "C", 0.0, 0.0, 0.0));
        }
        for i in 0..6 {
            mol.bonds
                .push(Bond::new(i, (i + 1) % 6, BondOrder::Aromatic));
        }

        let adj = AdjacencyList::from_molecule(&mol);
        for i in 0..6 {
            assert_eq!(adj.degree(i), 2);
            assert_eq!(adj.heavy_degree(i), 2);
        }
    }

    #[test]
    fn test_adjacency_list_empty() {
        let mol = Molecule::new("empty");
        let adj = AdjacencyList::from_molecule(&mol);
        assert_eq!(adj.num_atoms(), 0);
    }

    #[test]
    fn test_neighbor_atoms() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 2.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));

        let adj = AdjacencyList::from_molecule(&mol);
        let neighbors = adj.neighbor_atoms(1);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&0));
        assert!(neighbors.contains(&2));
    }

    #[test]
    fn test_bond_indices() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 2.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Double));

        let adj = AdjacencyList::from_molecule(&mol);
        let bonds = adj.bond_indices(1);
        assert_eq!(bonds.len(), 2);
        assert!(bonds.contains(&0));
        assert!(bonds.contains(&1));
    }

    #[test]
    fn test_is_hydrogen() {
        assert!(is_hydrogen("H"));
        assert!(is_hydrogen("D"));
        assert!(is_hydrogen("T"));
        assert!(is_hydrogen(" H "));
        assert!(!is_hydrogen("C"));
        assert!(!is_hydrogen("He"));
    }

    #[test]
    fn test_out_of_bounds() {
        let mol = Molecule::new("empty");
        let adj = AdjacencyList::from_molecule(&mol);
        assert_eq!(adj.degree(999), 0);
        assert_eq!(adj.heavy_degree(999), 0);
        assert!(adj.neighbors(999).is_empty());
    }
}
