//! OGB (Open Graph Benchmark) compatible GNN featurizer.
//!
//! Produces atom and bond feature matrices matching the OGB `AtomEncoder`/`BondEncoder`
//! conventions. A single call replaces the entire RDKit preprocessing pipeline.
//!
//! ## OGB Atom Features (9 features per atom)
//!
//! | Index | Feature | Encoding |
//! |-------|---------|----------|
//! | 0 | Atomic number | 1-118 |
//! | 1 | Chirality tag | 0=unspecified, 1=CW, 2=CCW, 3=other |
//! | 2 | Degree | 0-10 |
//! | 3 | Formal charge | shifted: charge + 5 (range 0-10) |
//! | 4 | Num Hs | 0-8 |
//! | 5 | Num radical electrons | 0-4 |
//! | 6 | Hybridization | 0=S, 1=SP, 2=SP2, 3=SP3, 4=SP3D, 5=SP3D2 |
//! | 7 | Is aromatic | 0 or 1 |
//! | 8 | Is in ring | 0 or 1 |
//!
//! ## OGB Bond Features (3 features per bond)
//!
//! | Index | Feature | Encoding |
//! |-------|---------|----------|
//! | 0 | Bond type | 0=single, 1=double, 2=triple, 3=aromatic |
//! | 1 | Bond stereo | 0=none, 1=any, 2=E/Z, 3=CIS, 4=TRANS |
//! | 2 | Is conjugated | 0 or 1 |
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//! use sdfrust::featurize::ogb;
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -0.3, 0.9, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! let atom_feats = ogb::ogb_atom_features(&mol);
//! assert_eq!(atom_feats.num_atoms, 3);
//! assert_eq!(atom_feats.features[0].len(), 9);
//!
//! let bond_feats = ogb::ogb_bond_features(&mol);
//! assert_eq!(bond_feats.num_bonds, 2);
//! ```

use crate::bond::{BondOrder, BondStereo};
use crate::descriptors::aromaticity::{all_aromatic_atoms, all_aromatic_bonds};
use crate::descriptors::chirality::all_chiralities;
use crate::descriptors::conjugation::all_conjugated_bonds;
use crate::descriptors::elements::get_element;
use crate::descriptors::hybridization::all_hybridizations;
use crate::descriptors::topological::ring_atoms;
use crate::descriptors::valence;
use crate::graph::AdjacencyList;
use crate::molecule::Molecule;

/// OGB atom feature matrix.
#[derive(Debug, Clone)]
pub struct OgbAtomFeatures {
    /// Number of atoms.
    pub num_atoms: usize,
    /// Feature matrix: `features[i]` is the 9-element feature vector for atom i.
    pub features: Vec<Vec<i32>>,
}

/// OGB bond feature matrix.
#[derive(Debug, Clone)]
pub struct OgbBondFeatures {
    /// Number of bonds.
    pub num_bonds: usize,
    /// Feature matrix: `features[i]` is the 3-element feature vector for bond i.
    pub features: Vec<Vec<i32>>,
}

/// Complete OGB graph features.
#[derive(Debug, Clone)]
pub struct OgbGraphFeatures {
    /// Atom features [N, 9].
    pub atom_features: OgbAtomFeatures,
    /// Bond features [E, 3] (directed, so E = 2 * num_bonds).
    pub bond_features: OgbBondFeatures,
    /// Edge index: (src, dst) arrays for directed edges.
    pub edge_src: Vec<usize>,
    /// Edge index: destination atoms.
    pub edge_dst: Vec<usize>,
}

/// Compute OGB-compatible atom features for a molecule.
///
/// Returns a matrix of shape `[N, 9]` where N is the number of atoms.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::featurize::ogb;
///
/// let mut mol = Molecule::new("methane");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
///
/// let feats = ogb::ogb_atom_features(&mol);
/// assert_eq!(feats.features[0][0], 6); // Carbon atomic number
/// assert_eq!(feats.features[1][0], 1); // Hydrogen atomic number
/// ```
pub fn ogb_atom_features(mol: &Molecule) -> OgbAtomFeatures {
    let n = mol.atom_count();
    let adj = AdjacencyList::from_molecule(mol);
    let hybridizations = all_hybridizations(mol);
    let aromatic_atoms = all_aromatic_atoms(mol);
    let in_ring = ring_atoms(mol);
    let chiralities = all_chiralities(mol);

    let mut features = Vec::with_capacity(n);

    for i in 0..n {
        let atom = &mol.atoms[i];
        let mut feat = vec![0i32; 9];

        // Feature 0: Atomic number
        feat[0] = get_element(&atom.element)
            .map(|e| e.atomic_number as i32)
            .unwrap_or(0);

        // Feature 1: Chirality tag (CIP-based perception)
        // 0=unspecified, 1=CW (R), 2=CCW (S), 3=other
        feat[1] = chiralities[i].to_ogb_index() as i32;

        // Feature 2: Degree
        feat[2] = adj.degree(i) as i32;

        // Feature 3: Formal charge (shifted: charge + 5)
        feat[3] = (atom.formal_charge as i32) + 5;

        // Feature 4: Number of Hs (implicit only, matching RDKit's GetTotalNumHs
        // when explicit H atoms are present as separate graph nodes)
        feat[4] = valence::implicit_hydrogen_count(mol, i) as i32;

        // Feature 5: Number of radical electrons
        feat[5] = atom.radical.unwrap_or(0) as i32;

        // Feature 6: Hybridization
        feat[6] = hybridizations[i].to_ogb_index() as i32;

        // Feature 7: Is aromatic
        feat[7] = if aromatic_atoms[i] { 1 } else { 0 };

        // Feature 8: Is in ring
        feat[8] = if in_ring[i] { 1 } else { 0 };

        features.push(feat);
    }

    OgbAtomFeatures {
        num_atoms: n,
        features,
    }
}

/// Compute OGB-compatible bond features for a molecule.
///
/// Returns a matrix of shape `[E, 3]` where E is the number of bonds.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::featurize::ogb;
///
/// let mut mol = Molecule::new("ethylene");
/// mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Double));
///
/// let feats = ogb::ogb_bond_features(&mol);
/// assert_eq!(feats.features[0][0], 1); // Double bond = 1
/// ```
pub fn ogb_bond_features(mol: &Molecule) -> OgbBondFeatures {
    let m = mol.bond_count();
    let conjugated = all_conjugated_bonds(mol);
    let aromatic_bonds = all_aromatic_bonds(mol);

    let mut features = Vec::with_capacity(m);

    for (i, bond) in mol.bonds.iter().enumerate() {
        let mut feat = vec![0i32; 3];

        // Feature 0: Bond type
        // 0=single, 1=double, 2=triple, 3=aromatic
        // Use perceived aromaticity to override Kekulized bond types
        feat[0] = if aromatic_bonds[i] {
            3 // Perceived aromatic
        } else {
            match bond.order {
                BondOrder::Single => 0,
                BondOrder::Double => 1,
                BondOrder::Triple => 2,
                BondOrder::Aromatic => 3,
                BondOrder::SingleOrDouble => 0,
                BondOrder::SingleOrAromatic => 0,
                BondOrder::DoubleOrAromatic => 1,
                BondOrder::Any => 0,
                BondOrder::Coordination => 0,
                BondOrder::Hydrogen => 0,
            }
        };

        // Feature 1: Bond stereo (E/Z double bond stereo only)
        // 0=none, 1=any, 2=E, 3=Z, 4=CIS, 5=TRANS
        // Note: SDF Up/Down on single bonds are chirality indicators (wedge bonds),
        // not E/Z stereo. Only encode stereo for double bonds.
        feat[1] = if bond.order == BondOrder::Double {
            match bond.stereo {
                BondStereo::None => 0,
                BondStereo::Up => 1, // E
                BondStereo::Either => 2,
                BondStereo::Down => 3, // Z
            }
        } else {
            0 // Single bonds: wedge stereo is for chirality, not OGB bond stereo
        };

        // Feature 2: Is conjugated
        feat[2] = if conjugated[i] { 1 } else { 0 };

        features.push(feat);
    }

    OgbBondFeatures {
        num_bonds: m,
        features,
    }
}

/// Compute complete OGB graph features including edge index.
///
/// Returns atom features, bond features, and a directed edge index
/// (both i→j and j→i for each bond) matching PyTorch Geometric's format.
///
/// # Example
///
/// ```rust
/// use sdfrust::{Molecule, Atom, Bond, BondOrder};
/// use sdfrust::featurize::ogb;
///
/// let mut mol = Molecule::new("water");
/// mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
/// mol.atoms.push(Atom::new(2, "H", -0.3, 0.9, 0.0));
/// mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
/// mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
///
/// let graph = ogb::ogb_graph_features(&mol);
/// assert_eq!(graph.atom_features.num_atoms, 3);
/// // Directed edges: 2 bonds × 2 = 4
/// assert_eq!(graph.edge_src.len(), 4);
/// ```
pub fn ogb_graph_features(mol: &Molecule) -> OgbGraphFeatures {
    let atom_features = ogb_atom_features(mol);
    let bond_feats = ogb_bond_features(mol);

    // Build directed edge index and duplicate bond features
    let mut edge_src = Vec::with_capacity(mol.bond_count() * 2);
    let mut edge_dst = Vec::with_capacity(mol.bond_count() * 2);
    let mut directed_bond_features = Vec::with_capacity(mol.bond_count() * 2);

    for (i, bond) in mol.bonds.iter().enumerate() {
        // Forward: atom1 → atom2
        edge_src.push(bond.atom1);
        edge_dst.push(bond.atom2);
        directed_bond_features.push(bond_feats.features[i].clone());

        // Reverse: atom2 → atom1
        edge_src.push(bond.atom2);
        edge_dst.push(bond.atom1);
        directed_bond_features.push(bond_feats.features[i].clone());
    }

    OgbGraphFeatures {
        atom_features,
        bond_features: OgbBondFeatures {
            num_bonds: directed_bond_features.len(),
            features: directed_bond_features,
        },
        edge_src,
        edge_dst,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::bond::{Bond, BondOrder};

    fn make_water() -> Molecule {
        let mut mol = Molecule::new("water");
        mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
        mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
        mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
        mol
    }

    fn make_benzene() -> Molecule {
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

    #[test]
    fn test_atom_features_water() {
        let mol = make_water();
        let feats = ogb_atom_features(&mol);

        assert_eq!(feats.num_atoms, 3);
        assert_eq!(feats.features.len(), 3);

        // O atom
        assert_eq!(feats.features[0][0], 8); // Atomic number of O
        assert_eq!(feats.features[0][2], 2); // Degree 2

        // H atoms
        assert_eq!(feats.features[1][0], 1); // Atomic number of H
        assert_eq!(feats.features[1][2], 1); // Degree 1
    }

    #[test]
    fn test_atom_features_9_elements() {
        let mol = make_water();
        let feats = ogb_atom_features(&mol);
        for feat in &feats.features {
            assert_eq!(feat.len(), 9);
        }
    }

    #[test]
    fn test_bond_features_water() {
        let mol = make_water();
        let feats = ogb_bond_features(&mol);

        assert_eq!(feats.num_bonds, 2);
        assert_eq!(feats.features.len(), 2);

        // Both are single bonds
        for feat in &feats.features {
            assert_eq!(feat[0], 0); // Single bond
            assert_eq!(feat.len(), 3);
        }
    }

    #[test]
    fn test_bond_features_benzene() {
        let mol = make_benzene();
        let feats = ogb_bond_features(&mol);

        for feat in &feats.features {
            assert_eq!(feat[0], 3); // Aromatic
            assert_eq!(feat[2], 1); // Conjugated
        }
    }

    #[test]
    fn test_atom_features_aromatic() {
        let mol = make_benzene();
        let feats = ogb_atom_features(&mol);

        for feat in &feats.features {
            assert_eq!(feat[7], 1); // Is aromatic
            assert_eq!(feat[8], 1); // Is in ring
        }
    }

    #[test]
    fn test_graph_features() {
        let mol = make_water();
        let graph = ogb_graph_features(&mol);

        assert_eq!(graph.atom_features.num_atoms, 3);
        // 2 bonds × 2 directions = 4 directed edges
        assert_eq!(graph.edge_src.len(), 4);
        assert_eq!(graph.edge_dst.len(), 4);
        assert_eq!(graph.bond_features.num_bonds, 4);
    }

    #[test]
    fn test_formal_charge_encoding() {
        let mut mol = Molecule::new("test");
        let mut atom = Atom::new(0, "N", 0.0, 0.0, 0.0);
        atom.formal_charge = 1;
        mol.atoms.push(atom);

        let feats = ogb_atom_features(&mol);
        assert_eq!(feats.features[0][3], 6); // 1 + 5 = 6
    }

    #[test]
    fn test_empty_molecule() {
        let mol = Molecule::new("empty");
        let atom_feats = ogb_atom_features(&mol);
        assert_eq!(atom_feats.num_atoms, 0);
        assert!(atom_feats.features.is_empty());

        let bond_feats = ogb_bond_features(&mol);
        assert_eq!(bond_feats.num_bonds, 0);
        assert!(bond_feats.features.is_empty());
    }

    #[test]
    fn test_hybridization_in_features() {
        // C with triple bond → SP → 1
        let mut mol = Molecule::new("test");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.2, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Triple));

        let feats = ogb_atom_features(&mol);
        assert_eq!(feats.features[0][6], 1); // SP = 1
    }
}
