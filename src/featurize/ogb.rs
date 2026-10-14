//! OGB (Open Graph Benchmark) compatible GNN featurizer.
//!
//! Produces atom and bond feature matrices matching the OGB `AtomEncoder`/`BondEncoder`
//! conventions. A single call replaces the entire RDKit preprocessing pipeline.
//!
//! Every feature is the index into the corresponding list of
//! `ogb.utils.features.allowable_features`, with out-of-range values mapped to the
//! trailing `misc` bucket (OGB's `safe_index`). The output can therefore be fed
//! directly into OGB's `AtomEncoder`/`BondEncoder`.
//!
//! ## OGB Atom Features (9 features per atom)
//!
//! | Index | Feature | Encoding (dim) |
//! |-------|---------|----------|
//! | 0 | Atomic number | Z - 1 (0..=117), misc = 118 (119) |
//! | 1 | Chirality tag | 0=unspecified, 1=CW, 2=CCW, 3=other (5) |
//! | 2 | Total degree (incl. implicit H) | 0..=10, misc = 11 (12) |
//! | 3 | Formal charge | charge + 5 for -5..=5, misc = 11 (12) |
//! | 4 | Total num Hs | 0..=8, misc = 9 (10) |
//! | 5 | Num radical electrons | 0..=4, misc = 5 (6) |
//! | 6 | Hybridization | 0=SP, 1=SP2, 2=SP3, 3=SP3D, 4=SP3D2, 5=misc (incl. S) (6) |
//! | 7 | Is aromatic | 0 or 1 (2) |
//! | 8 | Is in ring | 0 or 1 (2) |
//!
//! ## OGB Bond Features (3 features per bond)
//!
//! | Index | Feature | Encoding (dim) |
//! |-------|---------|----------|
//! | 0 | Bond type | 0=single, 1=double, 2=triple, 3=aromatic, 4=misc (5) |
//! | 1 | Bond stereo | 0=NONE, 1=Z, 2=E, 3=CIS, 4=TRANS, 5=ANY (6) |
//! | 2 | Is conjugated | 0 or 1 (2) |
//!
//! Bond stereo: E/Z perception from coordinates is not implemented yet, so only
//! `STEREONONE` (0) and `STEREOANY` (5, for "either" double bonds) are emitted.
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

/// OGB's `safe_index`: values outside `0..misc` map to the trailing `misc` bucket.
fn safe_index(value: i32, misc: i32) -> i32 {
    if (0..misc).contains(&value) {
        value
    } else {
        misc
    }
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
/// assert_eq!(feats.features[0][0], 5); // Carbon: Z - 1
/// assert_eq!(feats.features[1][0], 0); // Hydrogen: Z - 1
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

        // Feature 0: Atomic number, index into range(1, 119)
        feat[0] = get_element(&atom.element)
            .map(|e| safe_index(e.atomic_number as i32 - 1, 118))
            .unwrap_or(118);

        // Feature 1: Chirality tag (CHI_UNSPECIFIED, CW, CCW, OTHER)
        feat[1] = chiralities[i].to_ogb_index() as i32;

        let num_hs = valence::implicit_hydrogen_count(mol, i) as i32;

        // Feature 2: Total degree (explicit neighbors + implicit H), like GetTotalDegree
        feat[2] = safe_index(adj.degree(i) as i32 + num_hs, 11);

        // Feature 3: Formal charge, index into [-5..=5]
        feat[3] = safe_index(atom.formal_charge as i32 + 5, 11);

        // Feature 4: Total number of Hs on the atom (GetTotalNumHs; explicit H
        // atoms are separate graph nodes and counted in the degree instead)
        feat[4] = safe_index(num_hs, 9);

        // Feature 5: Number of radical electrons
        feat[5] = safe_index(atom.radical_electrons() as i32, 5);

        // Feature 6: Hybridization, index into [SP, SP2, SP3, SP3D, SP3D2]
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

        // Feature 0: Bond type, index into [SINGLE, DOUBLE, TRIPLE, AROMATIC]
        // Perceived aromaticity overrides Kekulé bond orders, as in RDKit.
        feat[0] = if aromatic_bonds[i] {
            3
        } else {
            match bond.order {
                BondOrder::Single => 0,
                BondOrder::Double => 1,
                BondOrder::Triple => 2,
                BondOrder::Aromatic => 3,
                _ => 4, // query, coordination and hydrogen bonds → misc
            }
        };

        // Feature 1: Bond stereo, index into
        // [STEREONONE, STEREOZ, STEREOE, STEREOCIS, STEREOTRANS, STEREOANY].
        // Wedge flags (Up/Down) describe tetrahedral centers, not E/Z, so only an
        // explicit "either" double bond is encoded (as STEREOANY).
        feat[1] = if bond.order == BondOrder::Double && bond.stereo == BondStereo::Either {
            5
        } else {
            0
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
        assert_eq!(feats.features[0][0], 7); // O: Z - 1
        assert_eq!(feats.features[0][2], 2); // Degree 2
        assert_eq!(feats.features[0][6], 2); // SP3

        // H atoms
        assert_eq!(feats.features[1][0], 0); // H: Z - 1
        assert_eq!(feats.features[1][2], 1); // Degree 1
        assert_eq!(feats.features[1][6], 5); // S hybridization → misc
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
        assert_eq!(feats.features[0][6], 0); // SP = 0
    }

    #[test]
    fn test_total_degree_includes_implicit_h() {
        // Heavy-atom-only methanol: C has 1 explicit neighbor + 3 implicit H
        let mut mol = Molecule::new("methanol");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "O", 1.4, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Single));

        let feats = ogb_atom_features(&mol);
        assert_eq!(feats.features[0][2], 4); // total degree
        assert_eq!(feats.features[0][4], 3); // num H
        assert_eq!(feats.features[1][2], 2);
        assert_eq!(feats.features[1][4], 1);
    }

    #[test]
    fn test_out_of_range_values_use_misc_bucket() {
        let mut mol = Molecule::new("odd");
        let mut a = Atom::new(0, "Xx", 0.0, 0.0, 0.0);
        a.formal_charge = 7;
        mol.atoms.push(a);
        mol.atoms.push(Atom::new(1, "C", 1.0, 0.0, 0.0));
        mol.bonds.push(Bond::new(0, 1, BondOrder::Coordination));

        let atoms = ogb_atom_features(&mol);
        assert_eq!(atoms.features[0][0], 118); // unknown element → misc
        assert_eq!(atoms.features[0][3], 11); // charge +7 → misc

        let bonds = ogb_bond_features(&mol);
        assert_eq!(bonds.features[0][0], 4); // coordination bond → misc
    }

    #[test]
    fn test_features_within_ogb_dims() {
        let atom_dims = [119, 5, 12, 12, 10, 6, 6, 2, 2];
        let bond_dims = [5, 6, 2];
        for mol in [make_water(), make_benzene()] {
            for f in ogb_atom_features(&mol).features {
                for (v, d) in f.iter().zip(atom_dims) {
                    assert!((0..d).contains(v));
                }
            }
            for f in ogb_bond_features(&mol).features {
                for (v, d) in f.iter().zip(bond_dims) {
                    assert!((0..d).contains(v));
                }
            }
        }
    }

    #[test]
    fn test_wedge_flags_on_double_bond_are_not_ez() {
        let mut mol = Molecule::new("ethylene");
        mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
        mol.atoms.push(Atom::new(1, "C", 1.3, 0.0, 0.0));
        let mut b = Bond::new(0, 1, BondOrder::Double);
        b.stereo = BondStereo::Up;
        mol.bonds.push(b);
        assert_eq!(ogb_bond_features(&mol).features[0][1], 0);

        mol.bonds[0].stereo = BondStereo::Either;
        assert_eq!(ogb_bond_features(&mol).features[0][1], 5); // STEREOANY
    }
}
