//! Extended Connectivity Fingerprints (ECFP) / Morgan fingerprints.
//!
//! Implements the Rogers & Hahn (2010) algorithm:
//! 1. Assign initial atom invariants (hash of atomic number, degree, etc.)
//! 2. Iteratively collect neighbor information at increasing radii
//! 3. Fold to a fixed-length bit vector
//!
//! This is a pure graph algorithm — no SMILES dependency required.
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
use crate::descriptors::elements::get_element;
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
/// * `n_bits` - Length of the folded bit vector (typically 1024 or 2048)
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

    // Fold to bit vector
    let mut bits = vec![false; n_bits];
    for &id in identifiers.values().flatten() {
        let bit_idx = (id as usize) % n_bits;
        bits[bit_idx] = true;
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
/// count-based similarity metrics or feature importance analysis.
///
/// # Arguments
///
/// * `mol` - The molecule
/// * `radius` - ECFP radius (2 for ECFP4, 3 for ECFP6)
pub fn ecfp_counts(mol: &Molecule, radius: usize) -> EcfpCountFingerprint {
    let identifiers = compute_ecfp_identifiers(mol, radius);

    let mut counts: HashMap<u32, u32> = HashMap::new();
    for &id in identifiers.values().flatten() {
        *counts.entry(id).or_insert(0) += 1;
    }

    EcfpCountFingerprint { counts, radius }
}

/// Core ECFP algorithm: compute identifiers for each atom at each radius.
fn compute_ecfp_identifiers(mol: &Molecule, radius: usize) -> HashMap<usize, Vec<u32>> {
    let n = mol.atom_count();
    if n == 0 {
        return HashMap::new();
    }

    let adj = AdjacencyList::from_molecule(mol);

    // Step 1: Initial atom invariants
    let mut identifiers: Vec<u32> = (0..n).map(|i| initial_atom_invariant(mol, i)).collect();

    // Collect all identifiers per atom
    let mut all_identifiers: HashMap<usize, Vec<u32>> = HashMap::new();
    for (i, &id) in identifiers.iter().enumerate() {
        all_identifiers.entry(i).or_default().push(id);
    }

    // Step 2: Iterative neighbor collection
    for _r in 0..radius {
        let mut new_identifiers = vec![0u32; n];

        for i in 0..n {
            let mut neighbor_info: Vec<(u32, u32)> = Vec::new();

            for &(neighbor, bond_idx) in adj.neighbors(i) {
                let bond_invariant = bond_type_hash(mol, bond_idx);
                neighbor_info.push((bond_invariant, identifiers[neighbor]));
            }

            // Sort for canonical ordering
            neighbor_info.sort();

            // Hash: combine current identifier with sorted neighbor info
            let mut hash = identifiers[i];
            for (bond_inv, atom_inv) in &neighbor_info {
                hash = hash_combine(hash, *bond_inv);
                hash = hash_combine(hash, *atom_inv);
            }

            new_identifiers[i] = hash;
            all_identifiers.entry(i).or_default().push(hash);
        }

        identifiers = new_identifiers;
    }

    all_identifiers
}

/// Compute initial atom invariant hash.
///
/// Based on Rogers & Hahn (2010): hash of atomic number, degree,
/// number of Hs, formal charge, ring membership.
fn initial_atom_invariant(mol: &Molecule, idx: usize) -> u32 {
    let atom = &mol.atoms[idx];
    let adj = AdjacencyList::from_molecule(mol);

    let atomic_num = get_element(&atom.element)
        .map(|e| e.atomic_number as u32)
        .unwrap_or(0);
    let degree = adj.degree(idx) as u32;
    let num_hs = valence::total_hydrogen_count(mol, idx) as u32;
    let charge = (atom.formal_charge as i32 + 5) as u32; // shift to positive
    let is_ring = if mol.is_atom_in_ring(idx) { 1u32 } else { 0 };

    let mut hash = atomic_num;
    hash = hash_combine(hash, degree);
    hash = hash_combine(hash, num_hs);
    hash = hash_combine(hash, charge);
    hash = hash_combine(hash, is_ring);

    hash
}

/// Hash for bond type.
fn bond_type_hash(mol: &Molecule, bond_idx: usize) -> u32 {
    match mol.bonds.get(bond_idx) {
        Some(bond) => match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
            _ => 0,
        },
        None => 0,
    }
}

/// Combine two hash values (FNV-1a inspired).
fn hash_combine(h: u32, val: u32) -> u32 {
    let mut hash = h;
    hash ^= val;
    hash = hash.wrapping_mul(16777619);
    hash
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
        assert!(sim >= 0.0 && sim <= 1.0);
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
