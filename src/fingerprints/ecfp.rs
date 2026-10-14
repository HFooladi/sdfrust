//! Extended Connectivity Fingerprints (ECFP) / Morgan fingerprints.
//!
//! Implements the Rogers & Hahn (2010) algorithm exactly as RDKit's Morgan
//! fingerprint generator does, so identifiers and folded bits are **identical
//! to RDKit** (`rdFingerprintGenerator.GetMorganGenerator(radius, fpSize)`
//! with default options) for the same molecular graph:
//! 1. Atom invariants: hash of atomic number, total degree, total H count,
//!    formal charge, isotope mass shift and ring membership
//! 2. Each layer hashes the layer index, the atom's identifier and its sorted
//!    (bond type, neighbour identifier) pairs; environments covering the same
//!    bonds as an earlier one are dropped
//! 3. Identifiers are folded into the bit vector modulo its length
//!
//! Aromaticity is perceived first, so Kekulé and aromatic inputs give the
//! same fingerprint. As in RDKit, explicit hydrogen atoms are part of the graph.
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::fingerprints::ecfp;
//!
//! let mut mol = Molecule::new("ethanol");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
//!
//! let fp = ecfp::ecfp(&mol, 2, 2048);
//! assert_eq!(fp.bits.len(), 2048);
//! assert!(fp.num_on_bits() > 0);
//! ```

use crate::bond::BondOrder;
use crate::descriptors::aromaticity::all_aromatic_bonds;
use crate::descriptors::elements::{element_by_atomic_number, get_element, isotope_mass};
use crate::descriptors::topological::ring_atoms;
use crate::descriptors::valence;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;
use std::collections::HashMap;

/// A folded bit-vector fingerprint.
#[derive(Debug, Clone)]
pub struct EcfpFingerprint {
    /// Bit vector (true = bit is set).
    pub bits: Vec<bool>,
    /// Number of bits in the fingerprint.
    pub n_bits: usize,
    /// Radius used for generation.
    pub radius: usize,
}

impl EcfpFingerprint {
    /// Returns the number of bits that are set (on).
    pub fn num_on_bits(&self) -> usize {
        self.bits.iter().filter(|&&b| b).count()
    }

    /// Returns the density (fraction of bits set).
    pub fn density(&self) -> f64 {
        if self.n_bits == 0 {
            return 0.0;
        }
        self.num_on_bits() as f64 / self.n_bits as f64
    }

    /// Returns the indices of set bits.
    pub fn on_bits(&self) -> Vec<usize> {
        self.bits
            .iter()
            .enumerate()
            .filter(|&(_, b)| *b)
            .map(|(i, _)| i)
            .collect()
    }

    /// Compute Tanimoto similarity to another fingerprint.
    ///
    /// Returns a value in [0, 1]. Returns 0.0 if both fingerprints are all zeros.
    pub fn tanimoto(&self, other: &EcfpFingerprint) -> f64 {
        let len = self.bits.len().min(other.bits.len());
        let mut a_count = 0u32;
        let mut b_count = 0u32;
        let mut ab_count = 0u32;

        for i in 0..len {
            if self.bits[i] {
                a_count += 1;
            }
            if other.bits[i] {
                b_count += 1;
            }
            if self.bits[i] && other.bits[i] {
                ab_count += 1;
            }
        }

        let denom = a_count + b_count - ab_count;
        if denom == 0 {
            return 0.0;
        }
        ab_count as f64 / denom as f64
    }
}

/// A count-based fingerprint (not folded to bits).
#[derive(Debug, Clone)]
pub struct EcfpCountFingerprint {
    /// Map from feature hash to count.
    pub counts: HashMap<u32, u32>,
    /// Radius used for generation.
    pub radius: usize,
}

/// Compute an ECFP fingerprint (bit vector).
///
/// # Arguments
///
/// * `mol` - The molecule
/// * `radius` - ECFP radius (2 for ECFP4, 3 for ECFP6)
/// * `n_bits` - Length of the folded bit vector (typically 1024 or 2048).
///   `0` yields an empty fingerprint.
///
/// Bits are identical to RDKit's Morgan fingerprint with default options
/// (`GetMorganGenerator(radius=radius, fpSize=n_bits).GetFingerprint(mol)`).
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::fingerprints::ecfp;
///
/// let mut mol = Molecule::new("methane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
///
/// let fp = ecfp::ecfp(&mol, 2, 1024);
/// assert_eq!(fp.n_bits, 1024);
/// ```
pub fn ecfp(mol: &Molecule, radius: usize, n_bits: usize) -> EcfpFingerprint {
    let identifiers = compute_ecfp_identifiers(mol, radius);

    // Fold to bit vector (n_bits == 0 yields an empty fingerprint)
    let mut bits = vec![false; n_bits];
    if n_bits > 0 {
        for &id in &identifiers {
            bits[(id as usize) % n_bits] = true;
        }
    }

    EcfpFingerprint {
        bits,
        n_bits,
        radius,
    }
}

/// Compute an ECFP count fingerprint (hash → count).
///
/// Returns unfolded feature identifiers with counts, useful for
/// count-based similarity metrics or feature importance analysis. Identical
/// to RDKit's `GetSparseCountFingerprint` from a default Morgan generator.
///
/// # Arguments
///
/// * `mol` - The molecule
/// * `radius` - ECFP radius (2 for ECFP4, 3 for ECFP6)
pub fn ecfp_counts(mol: &Molecule, radius: usize) -> EcfpCountFingerprint {
    let identifiers = compute_ecfp_identifiers(mol, radius);

    let mut counts: HashMap<u32, u32> = HashMap::new();
    for &id in &identifiers {
        *counts.entry(id).or_insert(0) += 1;
    }

    EcfpCountFingerprint { counts, radius }
}

/// Boost's `hash_combine` on a 32-bit seed, as used by RDKit's fingerprints.
fn hash_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// Boost's 64-bit `hash_combine` (used for `std::pair` and vector hashes).
fn hash_combine64(seed: u64, value: u64) -> u64 {
    seed ^ value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// RDKit bond type codes used as bond invariants.
fn bond_invariant(order: BondOrder, aromatic: bool) -> u32 {
    if aromatic {
        return 12;
    }
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 12,
        BondOrder::Hydrogen => 14,
        BondOrder::Coordination => 17,
        _ => 0,
    }
}

/// RDKit's connectivity invariants for every atom.
fn atom_invariants(mol: &Molecule, adj: &AdjacencyList) -> Vec<u32> {
    let implicit_h = valence::all_implicit_hydrogen_counts(mol);
    let in_ring = ring_atoms(mol);
    let atomic_number =
        |i: usize| get_element(&mol.atoms[i].element).map_or(0, |e| e.atomic_number);

    (0..mol.atom_count())
        .map(|i| {
            let atom = &mol.atoms[i];
            let element = get_element(&atom.element);
            let z = atomic_number(i);
            let h_neighbors = adj
                .neighbor_atoms(i)
                .into_iter()
                .filter(|&j| atomic_number(j) == 1)
                .count();
            let total_hs = implicit_h[i] as u32 + h_neighbors as u32;
            let total_degree = adj.degree(i) as u32 + implicit_h[i] as u32;
            // Mass shift for isotopes, truncated like C++'s static_cast<int>
            let delta_mass = element.map_or(0, |e| {
                let standard = element_by_atomic_number(e.atomic_number)
                    .map_or(e.atomic_weight, |s| s.atomic_weight);
                let mass = match atom.isotope {
                    Some(a) => isotope_mass(e.atomic_number, a).unwrap_or(a as f64),
                    None => e.atomic_weight, // D and T carry their isotope mass
                };
                (mass - standard) as i32
            });
            let mut components = vec![
                z as u32,
                total_degree,
                total_hs,
                atom.formal_charge as i32 as u32,
                delta_mass as u32,
            ];
            if in_ring[i] {
                components.push(1);
            }
            components.into_iter().fold(0u32, hash_combine)
        })
        .collect()
}

/// Core ECFP algorithm: the identifiers of every included environment
/// (layer 0 atoms plus each new, non-duplicate environment per layer).
fn compute_ecfp_identifiers(mol: &Molecule, radius: usize) -> Vec<u32> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    let adj = AdjacencyList::from_molecule(mol);
    let aromatic = all_aromatic_bonds(mol);
    let bond_inv: Vec<u32> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(b, bond)| bond_invariant(bond.order, aromatic[b]))
        .collect();

    let mut current = atom_invariants(mol, &adj);
    let mut identifiers = current.clone();

    // Bond sets covered by each atom's current environment
    let words = mol.bond_count().div_ceil(64).max(1);
    let mut neighborhoods = vec![vec![0u64; words]; n];
    let mut seen: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
    let mut dead = vec![false; n];

    for layer in 0..radius {
        let mut next = current.clone();
        let mut next_neighborhoods = neighborhoods.clone();
        let mut round: Vec<(Vec<u64>, u32, usize)> = Vec::new();

        for i in 0..n {
            if dead[i] {
                continue;
            }
            if adj.degree(i) == 0 {
                dead[i] = true;
                continue;
            }
            let mut env = neighborhoods[i].clone();
            let mut pairs: Vec<(u32, u32)> = Vec::with_capacity(adj.degree(i));
            for &(j, b) in adj.neighbors(i) {
                for (w, nw) in env.iter_mut().zip(&neighborhoods[j]) {
                    *w |= nw;
                }
                env[b / 64] |= 1 << (b % 64);
                pairs.push((bond_inv[b], current[j]));
            }
            pairs.sort_unstable();
            let mut id = hash_combine(layer as u32, current[i]);
            for (bt, nid) in pairs {
                // std::pair is hashed with Boost's 64-bit combine, then folded in
                let pair_hash = hash_combine64(hash_combine64(0, bt as u64), nid as u64);
                id = hash_combine(id, pair_hash as u32);
            }
            next[i] = id;
            next_neighborhoods[i] = env.clone();
            round.push((env, id, i));
        }

        // Equal environments: keep the smallest identifier, drop the rest
        round.sort();
        for (env, id, i) in round {
            if seen.contains(&env) {
                dead[i] = true;
            } else {
                identifiers.push(id);
                seen.insert(env);
            }
        }

        current = next;
        neighborhoods = next_neighborhoods;
    }

    identifiers
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    fn make_methane() -> Molecule {
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
        mol
    }

    fn make_ethanol() -> Molecule {
        let mut mol = Molecule::new("ethanol");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(1, 2, BondOrder::Single));
        mol
    }

    #[test]
    fn test_ecfp_basic() {
        let mol = make_methane();
        let fp = ecfp(&mol, 2, 2048);
        assert_eq!(fp.n_bits, 2048);
        assert!(fp.num_on_bits() > 0);
        assert!(fp.density() > 0.0);
        assert!(fp.density() < 1.0);
    }

    #[test]
    fn test_ecfp_radius_zero() {
        let mol = make_methane();
        let fp = ecfp(&mol, 0, 1024);
        assert!(fp.num_on_bits() > 0);
    }

    #[test]
    fn test_ecfp_same_molecule() {
        let mol1 = make_ethanol();
        let mol2 = make_ethanol();
        let fp1 = ecfp(&mol1, 2, 2048);
        let fp2 = ecfp(&mol2, 2, 2048);
        assert_eq!(fp1.tanimoto(&fp2), 1.0);
    }

    #[test]
    fn test_ecfp_different_molecules() {
        let mol1 = make_methane();
        let mol2 = make_ethanol();
        let fp1 = ecfp(&mol1, 2, 2048);
        let fp2 = ecfp(&mol2, 2, 2048);
        let sim = fp1.tanimoto(&fp2);
        assert!((0.0..=1.0).contains(&sim));
        // Different molecules should have < 1.0 similarity
        assert!(sim < 1.0);
    }

    #[test]
    fn test_ecfp_empty() {
        let mol = Molecule::new("empty");
        let fp = ecfp(&mol, 2, 1024);
        assert_eq!(fp.num_on_bits(), 0);
    }

    #[test]
    fn test_ecfp_counts() {
        let mol = make_ethanol();
        let fp = ecfp_counts(&mol, 2);
        assert!(!fp.counts.is_empty());
    }

    #[test]
    fn test_tanimoto_self() {
        let mol = make_ethanol();
        let fp = ecfp(&mol, 2, 2048);
        assert!((fp.tanimoto(&fp) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_tanimoto_empty() {
        let mol = Molecule::new("empty");
        let fp1 = ecfp(&mol, 2, 1024);
        let fp2 = ecfp(&mol, 2, 1024);
        assert_eq!(fp1.tanimoto(&fp2), 0.0);
    }

    #[test]
    fn test_zero_bits_does_not_panic() {
        let fp = ecfp(&make_ethanol(), 2, 0);
        assert_eq!(fp.n_bits, 0);
        assert_eq!(fp.num_on_bits(), 0);
    }

    #[test]
    fn test_on_bits() {
        let mol = make_methane();
        let fp = ecfp(&mol, 2, 2048);
        let on = fp.on_bits();
        assert_eq!(on.len(), fp.num_on_bits());
        for &idx in &on {
            assert!(fp.bits[idx]);
        }
    }

    #[test]
    fn test_ecfp_deterministic() {
        let mol = make_ethanol();
        let fp1 = ecfp(&mol, 2, 2048);
        let fp2 = ecfp(&mol, 2, 2048);
        assert_eq!(fp1.bits, fp2.bits);
    }

    #[test]
    fn test_higher_radius_more_bits() {
        let mol = make_ethanol();
        let fp0 = ecfp(&mol, 0, 2048);
        let fp2 = ecfp(&mol, 2, 2048);
        // Higher radius should generally set more (or equal) bits
        assert!(fp2.num_on_bits() >= fp0.num_on_bits());
    }
}
