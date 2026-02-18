//! ML-ready molecular featurization.
//!
//! Provides feature extraction functions for graph neural networks,
//! following OGB (Open Graph Benchmark) conventions.

pub mod ogb;

pub use ogb::{
    OgbAtomFeatures, OgbBondFeatures, ogb_atom_features, ogb_bond_features, ogb_graph_features,
};
