//! Cutoff-based neighbor list for 3D molecular models.
//!
//! Generates edge indices and pairwise distances in the format expected by
//! PyTorch Geometric and 3D equivariant GNNs (SchNet, DimeNet, PaiNN, GemNet).
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom};
//! use sdfrust::geometry::neighbor_list::{neighbor_list, NeighborList};
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
//!
//! let nl = neighbor_list(&mol, 2.0);
//! assert!(nl.num_edges() > 0);
//! ```

use crate::molecule::Molecule;

/// Result of a neighbor list computation.
///
/// Contains directed edges (both i→j and j→i) and pairwise distances,
/// matching the format expected by PyTorch Geometric's `Data` object.
#[derive(Debug, Clone)]
pub struct NeighborList {
    /// Source atom indices for each edge.
    pub edge_src: Vec<usize>,
    /// Destination atom indices for each edge.
    pub edge_dst: Vec<usize>,
    /// Euclidean distance for each edge.
    pub distances: Vec<f64>,
}

impl NeighborList {
    /// Returns the number of directed edges.
    pub fn num_edges(&self) -> usize {
        self.edge_src.len()
    }

    /// Returns the edge index as two separate vectors `(src, dst)`.
    ///
    /// This matches the PyTorch Geometric `edge_index` format when
    /// stacked as a `[2, E]` tensor.
    pub fn edge_index(&self) -> (&[usize], &[usize]) {
        (&self.edge_src, &self.edge_dst)
    }
}

/// Compute a neighbor list with a distance cutoff.
///
/// Returns directed pairs (both i→j and j→i) for all atom pairs within
/// the specified cutoff distance. This is the format required for GNN
/// message passing.
///
/// # Arguments
///
/// * `mol` - The molecule to compute the neighbor list for
/// * `cutoff` - Maximum distance in Angstroms
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::neighbor_list::neighbor_list;
///
/// let mut mol = Molecule::new("diatomic");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
///
/// let nl = neighbor_list(&mol, 2.0);
/// // Both i→j and j→i
/// assert_eq!(nl.num_edges(), 2);
/// assert!((nl.distances[0] - 1.5).abs() < 1e-10);
/// ```
pub fn neighbor_list(mol: &Molecule, cutoff: f64) -> NeighborList {
    let n = mol.atoms.len();
    let cutoff_sq = cutoff * cutoff;

    let mut edge_src = Vec::new();
    let mut edge_dst = Vec::new();
    let mut distances = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = mol.atoms[i].x - mol.atoms[j].x;
            let dy = mol.atoms[i].y - mol.atoms[j].y;
            let dz = mol.atoms[i].z - mol.atoms[j].z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq <= cutoff_sq {
                let dist = dist_sq.sqrt();
                // Add both directions
                edge_src.push(i);
                edge_dst.push(j);
                distances.push(dist);

                edge_src.push(j);
                edge_dst.push(i);
                distances.push(dist);
            }
        }
    }

    NeighborList {
        edge_src,
        edge_dst,
        distances,
    }
}

/// Compute a neighbor list with self-loops included.
///
/// Same as [`neighbor_list`] but also includes self-loop edges (i→i)
/// with distance 0.0. Some GNN architectures require self-loops.
///
/// # Arguments
///
/// * `mol` - The molecule to compute the neighbor list for
/// * `cutoff` - Maximum distance in Angstroms
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom};
/// use sdfrust::geometry::neighbor_list::neighbor_list_with_self_loops;
///
/// let mut mol = Molecule::new("single_atom");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
///
/// let nl = neighbor_list_with_self_loops(&mol, 5.0);
/// assert_eq!(nl.num_edges(), 1); // self-loop
/// assert!((nl.distances[0]).abs() < 1e-10);
/// ```
pub fn neighbor_list_with_self_loops(mol: &Molecule, cutoff: f64) -> NeighborList {
    let mut nl = neighbor_list(mol, cutoff);

    // Add self-loops
    for i in 0..mol.atoms.len() {
        nl.edge_src.push(i);
        nl.edge_dst.push(i);
        nl.distances.push(0.0);
    }

    nl
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;

    #[test]
    fn test_neighbor_list_two_atoms() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));

        let nl = neighbor_list(&mol, 2.0);
        assert_eq!(nl.num_edges(), 2);
        assert!((nl.distances[0] - 1.5).abs() < 1e-10);
        assert!((nl.distances[1] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_neighbor_list_cutoff_excludes() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 10.0, 0.0, 0.0));

        let nl = neighbor_list(&mol, 2.0);
        assert_eq!(nl.num_edges(), 0);
    }

    #[test]
    fn test_neighbor_list_triangle() {
        let mut mol = Molecule::new("triangle");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.5, 0.866, 0.0));

        // All pairs within ~1.0 Å
        let nl = neighbor_list(&mol, 1.5);
        assert_eq!(nl.num_edges(), 6); // 3 pairs × 2 directions
    }

    #[test]
    fn test_neighbor_list_empty() {
        let mol = Molecule::new("empty");
        let nl = neighbor_list(&mol, 5.0);
        assert_eq!(nl.num_edges(), 0);
    }

    #[test]
    fn test_neighbor_list_self_loops() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));

        let nl = neighbor_list_with_self_loops(&mol, 2.0);
        // 2 directed edges + 2 self-loops
        assert_eq!(nl.num_edges(), 4);
    }

    #[test]
    fn test_neighbor_list_edge_index() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.2, 0.0, 0.0));

        let nl = neighbor_list(&mol, 2.0);
        let (src, dst) = nl.edge_index();
        assert_eq!(src.len(), 2);
        assert_eq!(dst.len(), 2);
        // Should have 0→1 and 1→0
        assert!(src.contains(&0) && src.contains(&1));
    }

    #[test]
    fn test_neighbor_list_water() {
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.9572, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.2400, 0.9266, 0.0));

        // O-H distance ~ 0.96 Å, H-H distance ~ 1.51 Å
        let nl = neighbor_list(&mol, 1.0);
        assert_eq!(nl.num_edges(), 4); // O↔H1 and O↔H2

        let nl_wide = neighbor_list(&mol, 2.0);
        assert_eq!(nl_wide.num_edges(), 6); // all pairs
    }
}
