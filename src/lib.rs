//! # sdfrust
//!
//! A fast, pure-Rust parser for SDF (Structure Data File) and MOL2 chemical structure files.
//!
//! ## Features
//!
//! - Parse SDF V2000 format files (single and multi-molecule)
//! - Write SDF V2000 format files
//! - Iterate over large SDF files without loading everything into memory
//! - Access atom coordinates, bonds, and molecule properties
//! - Zero external dependencies for parsing (only `thiserror` for error handling)
//!
//! ## Quick Start
//!
//! ### Parse a single molecule
//!
//! ```rust,ignore
//! use sdfrust::{parse_sdf_file, Molecule};
//!
//! let molecule = parse_sdf_file("molecule.sdf")?;
//! println!("Name: {}", molecule.name);
//! println!("Atoms: {}", molecule.atom_count());
//! println!("Formula: {}", molecule.formula());
//! ```
//!
//! ### Parse multiple molecules
//!
//! ```rust,ignore
//! use sdfrust::parse_sdf_file_multi;
//!
//! let molecules = parse_sdf_file_multi("database.sdf")?;
//! for mol in &molecules {
//!     println!("{}: {} atoms", mol.name, mol.atom_count());
//! }
//! ```
//!
//! ### Iterate over a large file
//!
//! ```rust,ignore
//! use sdfrust::iter_sdf_file;
//!
//! for result in iter_sdf_file("large_database.sdf")? {
//!     let mol = result?;
//!     // Process each molecule without loading all into memory
//! }
//! ```
//!
//! ### Parse from string
//!
//! ```rust
//! use sdfrust::parse_sdf_string;
//!
//! let sdf_content = r#"methane
//!
//!
//!   5  4  0  0  0  0  0  0  0  0999 V2000
//!     0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//!     0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
//!    -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
//!    -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
//!     0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
//!   1  2  1  0  0  0  0
//!   1  3  1  0  0  0  0
//!   1  4  1  0  0  0  0
//!   1  5  1  0  0  0  0
//! M  END
//! $$$$
//! "#;
//!
//! let mol = parse_sdf_string(sdf_content).unwrap();
//! assert_eq!(mol.name, "methane");
//! assert_eq!(mol.atom_count(), 5);
//! assert_eq!(mol.formula(), "CH4");
//! ```
//!
//! ### Write molecules
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder, write_sdf_string};
//!
//! let mut mol = Molecule::new("water");
//! mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! let sdf_output = write_sdf_string(&mol).unwrap();
//! println!("{}", sdf_output);
//! ```

pub mod atom;
pub mod bond;
pub mod error;
pub mod molecule;
pub mod parser;
pub mod writer;

// Re-export main types
pub use atom::Atom;
pub use bond::{Bond, BondOrder, BondStereo};
pub use error::{Result, SdfError};
pub use molecule::Molecule;

// Re-export parser functions
pub use parser::{
    iter_sdf_file, parse_sdf_file, parse_sdf_file_multi, parse_sdf_string, parse_sdf_string_multi,
    SdfIterator, SdfParser,
};

// Re-export writer functions
pub use writer::{write_sdf, write_sdf_file, write_sdf_file_multi, write_sdf_multi, write_sdf_string};
