//! Molecular fingerprint algorithms.
//!
//! Provides fingerprint generation for similarity searching, ML features,
//! and virtual screening. Currently implements ECFP (Extended Connectivity
//! Fingerprints) / Morgan fingerprints.

pub mod ecfp;

pub use ecfp::{EcfpFingerprint, ecfp, ecfp_counts};
