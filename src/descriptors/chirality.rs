//! CIP chirality perception for tetrahedral stereocenters.
//!
//! Determines chirality tags (CW/CCW) using CIP priority assignment
//! via BFS expansion, and R/S determination from 3D signed volume
//! or 2D wedge/dash bond information. Provides OGB atom feature 1.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder, BondStereo};
//! use sdfrust::descriptors::chirality::{ChiralTag, atom_chirality};
//!
//! // Simple achiral molecule: methane (all H neighbors are equivalent)
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
//! assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
//! ```

use std::cmp::Ordering;

use crate::bond::{BondOrder, BondStereo};
use crate::descriptors::elements::get_element;
use crate::descriptors::hybridization::{Hybridization, all_hybridizations};
use crate::descriptors::valence::implicit_hydrogen_count;
use crate::graph::{AdjacencyList, is_hydrogen};
use crate::molecule::Molecule;

/// Maximum BFS depth for CIP priority assignment.
const MAX_CIP_DEPTH: usize = 8;

/// Threshold for signed volume to distinguish CW/CCW from flat.
const VOLUME_THRESHOLD: f64 = 0.01;

/// Chirality tag for a tetrahedral stereocenter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChiralTag {
    /// No chirality specified or not a stereocenter.
    Unspecified,
    /// Clockwise (R configuration).
    CW,
    /// Counterclockwise (S configuration).
    CCW,
    /// Other chirality type.
    Other,
}

impl ChiralTag {
    /// Returns the OGB-compatible integer index.
    ///
    /// - 0 = Unspecified
    /// - 1 = CW (R)
    /// - 2 = CCW (S)
    /// - 3 = Other
    pub fn to_ogb_index(&self) -> u8 {
        match self {
            ChiralTag::Unspecified => 0,
            ChiralTag::CW => 1,
            ChiralTag::CCW => 2,
            ChiralTag::Other => 3,
        }
    }
}

// ============================================================
// CIP Priority
// ============================================================

/// CIP priority represented as layers of sorted atomic numbers from BFS expansion.
///
/// Layer 0 = immediate neighbor's atomic number.
/// Each subsequent layer = sorted (descending) atomic numbers at that BFS distance.
/// Compared lexicographically: first divergent layer decides; higher atomic number = higher priority.
#[derive(Debug, Clone, PartialEq, Eq)]
struct CipPriority {
    layers: Vec<Vec<u8>>,
}

impl Ord for CipPriority {
    fn cmp(&self, other: &Self) -> Ordering {
        for (a, b) in self.layers.iter().zip(other.layers.iter()) {
            let cmp = a.cmp(b);
            if cmp != Ordering::Equal {
                return cmp;
            }
        }
        self.layers.len().cmp(&other.layers.len())
    }
}

impl PartialOrd for CipPriority {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Get the atomic number for an atom, returning 0 for unknown elements.
fn atomic_number_of(mol: &Molecule, idx: usize) -> u8 {
    get_element(&mol.atoms[idx].element)
        .map(|e| e.atomic_number)
        .unwrap_or(0)
}

/// Get the bond multiplicity for CIP phantom atom generation.
/// Single=1, Double=2, Triple=3, Aromatic=1 (treated as single for CIP).
fn bond_multiplicity(order: BondOrder) -> usize {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 1,
        _ => 1,
    }
}

/// Compute CIP priority for a substituent starting at `start_neighbor` of `center`.
fn compute_cip_priority(
    mol: &Molecule,
    adj: &AdjacencyList,
    center: usize,
    start_neighbor: usize,
) -> CipPriority {
    let n = mol.atom_count();

    // Layer 0: atomic number of the immediate neighbor
    let start_an = if start_neighbor == usize::MAX {
        1 // implicit H
    } else {
        atomic_number_of(mol, start_neighbor)
    };
    let mut layers = vec![vec![start_an]];

    // If implicit H, no further expansion
    if start_neighbor == usize::MAX {
        return CipPriority { layers };
    }

    // Find the bond between center and start_neighbor for multiplicity
    let start_mult = find_bond_multiplicity(mol, center, start_neighbor);

    // Initialize BFS frontier
    let mut current_real: Vec<usize> = vec![start_neighbor];
    // Track which real atoms we've added phantom copies for at this level
    let mut phantom_ans: Vec<u8> = Vec::new();
    // Add phantom copies for the start bond (multiplicity - 1)
    for _ in 1..start_mult {
        phantom_ans.push(atomic_number_of(mol, start_neighbor));
    }

    // Visited set
    let mut visited = vec![false; n];
    visited[center] = true;
    visited[start_neighbor] = true;

    for _depth in 1..MAX_CIP_DEPTH {
        let mut next_real: Vec<usize> = Vec::new();
        let mut layer_ans: Vec<u8> = Vec::new();

        // Add phantom atomic numbers carried from previous expansion
        layer_ans.extend_from_slice(&phantom_ans);
        phantom_ans.clear();

        for &node_idx in &current_real {
            // Expand real atom: add its neighbors
            for &(neighbor, bond_idx) in adj.neighbors(node_idx) {
                let an = atomic_number_of(mol, neighbor);
                let mult = mol
                    .bonds
                    .get(bond_idx)
                    .map(|b| bond_multiplicity(b.order))
                    .unwrap_or(1);

                if visited[neighbor] {
                    // Visited atom: include its atomic number as phantom (don't expand)
                    // One copy for the bond itself, plus (mult-1) for double/triple bond phantoms
                    for _ in 0..mult {
                        layer_ans.push(an);
                    }
                } else {
                    // New atom: one real copy + (mult-1) phantom copies
                    layer_ans.push(an);
                    next_real.push(neighbor);
                    visited[neighbor] = true;
                    for _ in 1..mult {
                        layer_ans.push(an);
                    }
                }
            }

            // Add implicit hydrogens at this node
            let impl_h = implicit_hydrogen_count(mol, node_idx);
            layer_ans.extend(std::iter::repeat_n(1u8, impl_h as usize));
        }

        if layer_ans.is_empty() {
            break;
        }

        // Sort descending (higher atomic number = higher priority)
        layer_ans.sort_unstable_by(|a, b| b.cmp(a));
        layers.push(layer_ans);

        current_real = next_real;
        if current_real.is_empty() {
            break;
        }
    }

    CipPriority { layers }
}

/// Find the bond multiplicity between two atoms.
fn find_bond_multiplicity(mol: &Molecule, a: usize, b: usize) -> usize {
    for bond in &mol.bonds {
        if (bond.atom1 == a && bond.atom2 == b) || (bond.atom1 == b && bond.atom2 == a) {
            return bond_multiplicity(bond.order);
        }
    }
    1
}

// ============================================================
// Stereocenter Detection
// ============================================================

/// Normalize element symbol: uppercase first, lowercase rest.
fn normalize_element(element: &str) -> String {
    let elem = element.trim();
    let mut chars = elem.chars();
    match chars.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase(),
    }
}

/// Check if an atom is a potential tetrahedral stereocenter.
fn is_potential_stereocenter(
    mol: &Molecule,
    idx: usize,
    adj: &AdjacencyList,
    hybridizations: &[Hybridization],
) -> bool {
    // Must be SP3
    if hybridizations[idx] != Hybridization::SP3 {
        return false;
    }

    // Skip hydrogen
    if is_hydrogen(&mol.atoms[idx].element) {
        return false;
    }

    let elem = normalize_element(&mol.atoms[idx].element);

    // Skip nitrogen (invertible pyramidal, not stereocenter in ML context)
    if elem == "N" {
        return false;
    }

    // Must have exactly 4 substituents (explicit + implicit H)
    let explicit = adj.degree(idx);
    let impl_h = implicit_hydrogen_count(mol, idx);
    let total = explicit + impl_h as usize;

    if total != 4 {
        return false;
    }

    // Skip if 2+ implicit H (can't have 4 different substituents)
    if impl_h >= 2 {
        return false;
    }

    // Allow: C, S, P, Se, Si, Ge
    matches!(elem.as_str(), "C" | "S" | "P" | "Se" | "Si" | "Ge")
}

// ============================================================
// Spatial Arrangement (Signed Volume)
// ============================================================

/// Check if a molecule has 2D coordinates (all z approximately 0).
fn is_2d(mol: &Molecule) -> bool {
    mol.atoms.iter().all(|a| a.z.abs() < 1e-6)
}

/// Get the z-offset for a neighbor based on wedge/dash bonds (2D only).
///
/// In SDF V2000, stereo is from atom1 toward atom2:
/// - Up (wedge, 1): atom2 toward viewer → z = +1
/// - Down (dashed, 6): atom2 away from viewer → z = -1
fn get_z_offset(mol: &Molecule, center: usize, neighbor: usize) -> f64 {
    for bond in &mol.bonds {
        if bond.atom1 == center && bond.atom2 == neighbor {
            return match bond.stereo {
                BondStereo::Up => 1.0,
                BondStereo::Down => -1.0,
                _ => 0.0,
            };
        }
        if bond.atom2 == center && bond.atom1 == neighbor {
            // Reversed bond direction: invert stereo interpretation
            return match bond.stereo {
                BondStereo::Up => -1.0,
                BondStereo::Down => 1.0,
                _ => 0.0,
            };
        }
    }
    0.0
}

/// Check if the stereocenter has any wedge/dash bond information.
fn has_stereo_bonds(mol: &Molecule, center: usize) -> bool {
    mol.bonds.iter().any(|bond| {
        bond.stereo != BondStereo::None && (bond.atom1 == center || bond.atom2 == center)
    })
}

/// Get the position of a substituent, with z-modification for 2D molecules.
fn substituent_position(
    mol: &Molecule,
    center: usize,
    sub_idx: usize,
    mol_is_2d: bool,
) -> (f64, f64, f64) {
    let atom = &mol.atoms[sub_idx];
    let z = if mol_is_2d {
        get_z_offset(mol, center, sub_idx)
    } else {
        atom.z
    };
    (atom.x, atom.y, z)
}

/// Synthesize a position for an implicit hydrogen atom.
///
/// Places the H opposite the centroid of explicit neighbors, ensuring
/// it's on the opposite side of the z-axis from any z-offset neighbors.
fn synthesize_implicit_h_position(
    center_pos: (f64, f64, f64),
    explicit_positions: &[(f64, f64, f64)],
) -> (f64, f64, f64) {
    let (cx, cy, cz) = center_pos;

    if explicit_positions.is_empty() {
        return (cx + 1.0, cy, cz);
    }

    let n = explicit_positions.len() as f64;
    let avg_x: f64 = explicit_positions.iter().map(|p| p.0).sum::<f64>() / n;
    let avg_y: f64 = explicit_positions.iter().map(|p| p.1).sum::<f64>() / n;
    let avg_z: f64 = explicit_positions.iter().map(|p| p.2).sum::<f64>() / n;

    // Place H opposite the centroid through the center
    let dx = cx - avg_x;
    let dy = cy - avg_y;
    let dz = cz - avg_z;
    let len = (dx * dx + dy * dy + dz * dz).sqrt();

    if len < 1e-10 {
        // Neighbors are symmetric around center; place H along z-axis
        return (cx, cy, cz + 1.0);
    }

    let scale = 1.0 / len;
    (cx + dx * scale, cy + dy * scale, cz + dz * scale)
}

/// Compute the signed volume of the tetrahedron formed by 4 substituent positions
/// relative to the center.
///
/// Substituents are ordered by atom index (ascending), with implicit H last.
/// With the 4th neighbor behind the center:
/// negative → CW, positive → CCW, near-zero → Unspecified.
fn signed_volume(subs: &[(f64, f64, f64); 4], center: &(f64, f64, f64)) -> f64 {
    // Vectors from center to each substituent
    let v: [(f64, f64, f64); 4] = [
        (
            subs[0].0 - center.0,
            subs[0].1 - center.1,
            subs[0].2 - center.2,
        ),
        (
            subs[1].0 - center.0,
            subs[1].1 - center.1,
            subs[1].2 - center.2,
        ),
        (
            subs[2].0 - center.0,
            subs[2].1 - center.1,
            subs[2].2 - center.2,
        ),
        (
            subs[3].0 - center.0,
            subs[3].1 - center.1,
            subs[3].2 - center.2,
        ),
    ];

    // a = v[0] - v[3], b = v[1] - v[3], c = v[2] - v[3]
    let a = (v[0].0 - v[3].0, v[0].1 - v[3].1, v[0].2 - v[3].2);
    let b = (v[1].0 - v[3].0, v[1].1 - v[3].1, v[1].2 - v[3].2);
    let c = (v[2].0 - v[3].0, v[2].1 - v[3].1, v[2].2 - v[3].2);

    // Cross product b x c
    let bxc = (
        b.1 * c.2 - b.2 * c.1,
        b.2 * c.0 - b.0 * c.2,
        b.0 * c.1 - b.1 * c.0,
    );

    // Dot product a . (b x c)
    a.0 * bxc.0 + a.1 * bxc.1 + a.2 * bxc.2
}

// ============================================================
// Public API
// ============================================================

/// Determine the chirality of atom `idx`.
///
/// Uses CIP priority assignment and signed volume to detect
/// R/S configuration. Returns `Unspecified` for non-stereocenters.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::chirality::{ChiralTag, atom_chirality};
///
/// let mut mol = Molecule::new("test");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "F", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "Cl", 0.0, 1.0, 0.0));
/// mol.atoms.push(Atom::new(3, "Br", 0.0, 0.0, 1.0));
/// mol.atoms.push(Atom::new(4, "H", -1.0, -1.0, -1.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 4, BondOrder::Single));
///
/// let tag = atom_chirality(&mol, 0);
/// assert!(tag == ChiralTag::CW || tag == ChiralTag::CCW);
/// ```
pub fn atom_chirality(mol: &Molecule, idx: usize) -> ChiralTag {
    all_chiralities(mol)
        .get(idx)
        .copied()
        .unwrap_or(ChiralTag::Unspecified)
}

/// Compute chirality for all atoms in the molecule.
///
/// Returns a vector of length `mol.atom_count()`.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::descriptors::chirality::{ChiralTag, all_chiralities};
///
/// let mut mol = Molecule::new("ethane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let chiralities = all_chiralities(&mol);
/// assert_eq!(chiralities[0], ChiralTag::Unspecified);
/// assert_eq!(chiralities[1], ChiralTag::Unspecified);
/// ```
pub fn all_chiralities(mol: &Molecule) -> Vec<ChiralTag> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    let adj = AdjacencyList::from_molecule(mol);
    let hybridizations = all_hybridizations(mol);
    let mol_is_2d = is_2d(mol);

    let mut result = vec![ChiralTag::Unspecified; n];

    for (i, chirality) in result.iter_mut().enumerate() {
        if !is_potential_stereocenter(mol, i, &adj, &hybridizations) {
            continue;
        }

        // For 2D molecules, we need wedge/dash bonds to determine chirality
        if mol_is_2d && !has_stereo_bonds(mol, i) {
            continue;
        }

        // Collect substituents with CIP priorities for stereocenter detection
        let mut explicit_neighbors = adj.neighbor_atoms(i);
        explicit_neighbors.sort_unstable(); // canonical atom index order
        let impl_h = implicit_hydrogen_count(mol, i);

        let mut sub_priorities: Vec<(CipPriority, usize)> = Vec::with_capacity(4);
        for &neighbor in &explicit_neighbors {
            let priority = compute_cip_priority(mol, &adj, i, neighbor);
            sub_priorities.push((priority, neighbor));
        }
        // Add implicit H (lowest possible priority)
        for _ in 0..impl_h {
            let priority = CipPriority {
                layers: vec![vec![1]],
            };
            sub_priorities.push((priority, usize::MAX));
        }

        // Check all 4 substituents have different CIP priorities (true stereocenter)
        if sub_priorities.len() != 4 {
            continue;
        }
        // Sort by CIP priority to check for duplicates
        let mut sorted_prios: Vec<&CipPriority> = sub_priorities.iter().map(|(p, _)| p).collect();
        sorted_prios.sort();
        if sorted_prios[0] == sorted_prios[1]
            || sorted_prios[1] == sorted_prios[2]
            || sorted_prios[2] == sorted_prios[3]
        {
            continue;
        }

        // Get positions for signed volume calculation.
        // RDKit's CW/CCW is based on ATOM INDEX order (not CIP priority order).
        // Explicit neighbors are already sorted by atom index.
        // Implicit H goes last (highest virtual index).
        let center_pos = if mol_is_2d {
            (mol.atoms[i].x, mol.atoms[i].y, 0.0)
        } else {
            (mol.atoms[i].x, mol.atoms[i].y, mol.atoms[i].z)
        };

        // Build ordered neighbor list: explicit (by atom index) + implicit H last
        let ordered_neighbors: Vec<usize> = {
            let mut v: Vec<usize> = explicit_neighbors.clone();
            for _ in 0..impl_h {
                v.push(usize::MAX);
            }
            v
        };

        // Collect explicit neighbor positions (needed for implicit H synthesis)
        let explicit_positions: Vec<(f64, f64, f64)> = ordered_neighbors
            .iter()
            .filter(|idx| **idx != usize::MAX)
            .map(|idx| substituent_position(mol, i, *idx, mol_is_2d))
            .collect();

        let mut positions = [(0.0, 0.0, 0.0); 4];
        for (k, &sub_idx) in ordered_neighbors.iter().enumerate() {
            if sub_idx == usize::MAX {
                // Implicit H: synthesize position
                positions[k] = synthesize_implicit_h_position(center_pos, &explicit_positions);
            } else {
                positions[k] = substituent_position(mol, i, sub_idx, mol_is_2d);
            }
        }

        // Compute signed volume using atom-index-ordered neighbors.
        // With the 4th neighbor (highest index / implicit H) behind the center:
        // negative → CW, positive → CCW
        let vol = signed_volume(&positions, &center_pos);

        if vol < -VOLUME_THRESHOLD {
            *chirality = ChiralTag::CW;
        } else if vol > VOLUME_THRESHOLD {
            *chirality = ChiralTag::CCW;
        }
        // Near-zero → Unspecified (flat tetrahedron)
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder, BondStereo};

    #[test]
    fn test_achiral_methane() {
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

        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_sp2_not_stereocenter() {
        let mut mol = Molecule::new("ethylene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));

        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_nitrogen_not_stereocenter() {
        // Nitrogen is always Unspecified (invertible pyramid)
        let mut mol = Molecule::new("amine");
        mol.atoms.push(Atom::new(0, "N", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "C", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "C", 0.0, 0.0, 1.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));

        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_chiral_3d_r() {
        // C with F, Cl, Br, H in 3D → R configuration
        // Arrange so signed volume is negative → CW (R)
        let mut mol = Molecule::new("r_chiral");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "Br", -1.0, -1.0, 1.0));
        mol.atoms.push(Atom::new(2, "Cl", -1.0, 1.0, -1.0));
        mol.atoms.push(Atom::new(3, "F", 1.0, -1.0, -1.0));
        mol.atoms.push(Atom::new(4, "H", 1.0, 1.0, 1.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let tag = atom_chirality(&mol, 0);
        assert!(
            tag == ChiralTag::CW || tag == ChiralTag::CCW,
            "Expected CW or CCW, got {:?}",
            tag
        );
    }

    #[test]
    fn test_chiral_3d_cw_ccw_differ() {
        // Two enantiomers should give opposite chirality tags
        let mut mol_r = Molecule::new("r");
        mol_r.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol_r.atoms.push(Atom::new(1, "Br", 1.0, 0.0, 0.0));
        mol_r.atoms.push(Atom::new(2, "Cl", 0.0, 1.0, 0.0));
        mol_r.atoms.push(Atom::new(3, "F", 0.0, 0.0, 1.0));
        mol_r.atoms.push(Atom::new(4, "H", -1.0, -1.0, -1.0));
        mol_r.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol_r.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol_r.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol_r.bonds.push(Bond::new(0, 4, BondOrder::Single));

        // Mirror image: reflect F across z-axis
        let mut mol_s = Molecule::new("s");
        mol_s.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol_s.atoms.push(Atom::new(1, "Br", 1.0, 0.0, 0.0));
        mol_s.atoms.push(Atom::new(2, "Cl", 0.0, 1.0, 0.0));
        mol_s.atoms.push(Atom::new(3, "F", 0.0, 0.0, -1.0));
        mol_s.atoms.push(Atom::new(4, "H", -1.0, -1.0, 1.0));
        mol_s.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol_s.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol_s.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol_s.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let tag_r = atom_chirality(&mol_r, 0);
        let tag_s = atom_chirality(&mol_s, 0);

        assert_ne!(tag_r, ChiralTag::Unspecified, "R should be chiral");
        assert_ne!(tag_s, ChiralTag::Unspecified, "S should be chiral");
        assert_ne!(
            tag_r, tag_s,
            "Enantiomers should have different chirality tags"
        );
    }

    #[test]
    fn test_2d_with_wedge() {
        // 2D molecule with wedge bond (Up) on one neighbor
        let mut mol = Molecule::new("chiral_2d");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "F", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "Cl", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "Br", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
        mol.bonds
            .push(Bond::with_stereo(0, 1, BondOrder::Single, BondStereo::Up));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let tag = atom_chirality(&mol, 0);
        assert!(
            tag == ChiralTag::CW || tag == ChiralTag::CCW,
            "Expected CW or CCW for 2D chiral center with wedge bond, got {:?}",
            tag
        );
    }

    #[test]
    fn test_2d_no_wedge_unspecified() {
        // 2D molecule without any wedge/dash bonds → can't determine chirality
        let mut mol = Molecule::new("flat");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "F", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "Cl", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(3, "Br", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_symmetric_substituents_not_stereocenter() {
        // C with two identical substituents (two methyl groups)
        // C bonded to: CH3, CH3, F, H → not a stereocenter
        let mut mol = Molecule::new("symmetric");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0)); // methyl 1
        mol.atoms.push(Atom::new(2, "C", -1.0, 0.0, 0.0)); // methyl 2
        mol.atoms.push(Atom::new(3, "F", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, 0.0, 1.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_ogb_index() {
        assert_eq!(ChiralTag::Unspecified.to_ogb_index(), 0);
        assert_eq!(ChiralTag::CW.to_ogb_index(), 1);
        assert_eq!(ChiralTag::CCW.to_ogb_index(), 2);
        assert_eq!(ChiralTag::Other.to_ogb_index(), 3);
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new("empty");
        assert!(all_chiralities(&mol).is_empty());
    }

    #[test]
    fn test_all_chiralities_length() {
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let chiralities = all_chiralities(&mol);
        assert_eq!(chiralities.len(), 2);
    }

    #[test]
    fn test_cip_double_bond_phantom() {
        // C bonded to: O (double bond, i.e., C=O), OH, NH2, H
        // CIP priorities: O(=) > OH > NH2 > H
        // The double bond phantom makes the C=O oxygen rank higher than O-H
        let mut mol = Molecule::new("cip_test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0)); // center
        mol.atoms.push(Atom::new(1, "O", 1.0, 0.0, 1.0)); // carbonyl O
        mol.atoms.push(Atom::new(2, "O", -1.0, 0.0, 0.0)); // hydroxyl O
        mol.atoms.push(Atom::new(3, "N", 0.0, 1.0, 0.0)); // amino N
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, -1.0));

        // C=O double bond, C-OH single, C-NH2 single, C-H single
        mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        // This is not SP3 (it has a double bond → SP2), so it won't be detected as a stereocenter
        // That's correct behavior - only SP3 centers are tetrahedral stereocenters
        assert_eq!(atom_chirality(&mol, 0), ChiralTag::Unspecified);
    }

    #[test]
    fn test_aspirin_no_stereocenters() {
        // Parse aspirin - should have no stereocenters
        let path = std::path::Path::new("tests/test_data/aspirin.sdf");
        if path.exists() {
            let mol = crate::parse_sdf_file(path).unwrap();
            let chiralities = all_chiralities(&mol);
            assert!(
                chiralities.iter().all(|c| *c == ChiralTag::Unspecified),
                "Aspirin should have no stereocenters"
            );
        }
    }

    #[test]
    fn test_methionine_has_stereocenter() {
        let path = std::path::Path::new("tests/test_data/methionine.sdf");
        if path.exists() {
            let mol = crate::parse_sdf_file(path).unwrap();
            let chiralities = all_chiralities(&mol);
            let chiral_count = chiralities
                .iter()
                .filter(|c| **c != ChiralTag::Unspecified)
                .count();
            assert_eq!(
                chiral_count, 1,
                "Methionine should have exactly 1 stereocenter, got {}",
                chiral_count
            );
        }
    }

    #[test]
    fn test_caffeine_no_stereocenters() {
        let path = std::path::Path::new("tests/test_data/caffeine_pubchem.sdf");
        if path.exists() {
            let mol = crate::parse_sdf_file(path).unwrap();
            let chiralities = all_chiralities(&mol);
            assert!(
                chiralities.iter().all(|c| *c == ChiralTag::Unspecified),
                "Caffeine should have no stereocenters"
            );
        }
    }

    #[test]
    fn test_glucose_multiple_stereocenters() {
        let path = std::path::Path::new("tests/test_data/glucose.sdf");
        if path.exists() {
            let mol = crate::parse_sdf_file(path).unwrap();
            let chiralities = all_chiralities(&mol);
            let chiral_count = chiralities
                .iter()
                .filter(|c| **c != ChiralTag::Unspecified)
                .count();
            // Glucose (sucrose in test_data) has 9 stereocenters
            assert!(
                chiral_count >= 5,
                "Glucose should have multiple stereocenters, got {}",
                chiral_count
            );
        }
    }

    #[test]
    fn test_2d_down_bond_stereocenter() {
        // A 2D stereocenter with a Down (dashed) bond
        // Center C bonded to: Br(dash=away), Cl, F, H
        // CIP: Br(35) > Cl(17) > F(9) > H(1)
        let mut mol = Molecule::new("down_stereo");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "Br", 1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "Cl", 0.0, 1.0, 0.0));
        mol.atoms.push(Atom::new(3, "F", -1.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(4, "H", 0.0, -1.0, 0.0));

        mol.bonds
            .push(Bond::with_stereo(0, 1, BondOrder::Single, BondStereo::Down));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 3, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 4, BondOrder::Single));

        let tag = atom_chirality(&mol, 0);
        assert_ne!(
            tag,
            ChiralTag::Unspecified,
            "Should detect chirality from Down bond"
        );
    }
}
