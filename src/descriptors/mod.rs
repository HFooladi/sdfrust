//! Molecular descriptors for chemical structure analysis.
//!
//! This module provides functions for calculating molecular properties:
//!
//! - **Molecular properties**: molecular weight, exact mass, heavy atom count
//! - **Topological properties**: ring count, rotatable bonds
//!
//! # Example
//!
//! ```rust
//! use sdfrust::{parse_sdf_string, descriptors};
//!
//! let sdf = r#"water
//!
//!
//!   3  2  0  0  0  0  0  0  0  0999 V2000
//!     0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
//!     0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
//!    -0.2400    0.9266    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
//!   1  2  1  0  0  0  0
//!   1  3  1  0  0  0  0
//! M  END
//! $$$$
//! "#;
//!
//! let mol = parse_sdf_string(sdf).unwrap();
//! let mw = descriptors::molecular_weight(&mol).unwrap();
//! assert!((mw - 18.015).abs() < 0.01);
//! ```

pub mod elements;
pub mod molecular;
pub mod topological;

// Re-export element data functions
pub use elements::{atomic_weight, get_element, monoisotopic_mass, ElementData};

// Re-export molecular descriptor functions
pub use molecular::{bond_type_counts, exact_mass, heavy_atom_count, molecular_weight};

// Re-export topological descriptor functions
pub use topological::{ring_atoms, ring_bonds, ring_count, rotatable_bond_count};
