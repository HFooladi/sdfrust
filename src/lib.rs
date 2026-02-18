//! # sdfrust
//!
//! A fast, pure-Rust parser for SDF (Structure Data File), MOL2, and XYZ chemical structure files.
//!
//! ## Features
//!
//! - Parse SDF V2000 and V3000 format files (single and multi-molecule)
//! - Parse TRIPOS MOL2 format files (single and multi-molecule)
//! - Parse XYZ coordinate files (single and multi-molecule, atomic numbers supported)
//! - Write SDF V2000 and V3000 format files
//! - Automatic format detection for SDF, MOL2, and XYZ files
//! - Transparent gzip decompression for all file parsers (optional `gzip` feature)
//! - Support for molecules with >999 atoms/bonds (V3000)
//! - Enhanced stereochemistry, SGroups, and collections (V3000)
//! - Iterate over large files without loading everything into memory
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
//!
//! ### Parse MOL2 files
//!
//! ```rust
//! use sdfrust::parse_mol2_string;
//!
//! let mol2_content = r#"@<TRIPOS>MOLECULE
//! water
//!  3 2 0 0 0
//! SMALL
//! NO_CHARGES
//!
//! @<TRIPOS>ATOM
//!       1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
//!       2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
//!       3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
//! @<TRIPOS>BOND
//!      1     1     2 1
//!      2     1     3 1
//! "#;
//!
//! let mol = parse_mol2_string(mol2_content).unwrap();
//! assert_eq!(mol.name, "water");
//! assert_eq!(mol.formula(), "H2O");
//! ```
//!
//! ## Error Handling
//!
//! All parsing functions return `Result<T, SdfError>`. The library provides specific
//! error variants for different failure modes:
//!
//! ```rust
//! use sdfrust::{parse_sdf_string, SdfError};
//!
//! let result = parse_sdf_string("invalid content");
//! match result {
//!     Ok(mol) => println!("Parsed: {}", mol.name),
//!     Err(SdfError::EmptyFile) => println!("File was empty"),
//!     Err(SdfError::Parse { line, message }) => {
//!         println!("Parse error at line {}: {}", line, message);
//!     }
//!     Err(SdfError::InvalidCountsLine(s)) => {
//!         println!("Bad counts line: {}", s);
//!     }
//!     Err(e) => println!("Other error: {}", e),
//! }
//! ```
//!
//! ### Common Error Types
//!
//! - `SdfError::Io` - File I/O errors (file not found, permission denied)
//! - `SdfError::Parse` - General parse errors with line number
//! - `SdfError::EmptyFile` - The file contains no data
//! - `SdfError::AtomCountMismatch` - Declared atom count doesn't match actual atoms
//! - `SdfError::BondCountMismatch` - Declared bond count doesn't match actual bonds
//! - `SdfError::InvalidAtomIndex` - Bond references non-existent atom
//! - `SdfError::InvalidBondOrder` - Unrecognized bond type
//! - `SdfError::InvalidCountsLine` - Malformed counts line in header
//! - `SdfError::MissingSection` - Required section not found (MOL2)
//!
//! ### Handling Multi-Molecule Files with Errors
//!
//! When iterating, each molecule is parsed independently:
//!
//! ```rust,ignore
//! use sdfrust::iter_sdf_file;
//!
//! let mut success_count = 0;
//! let mut error_count = 0;
//!
//! for result in iter_sdf_file("database.sdf")? {
//!     match result {
//!         Ok(mol) => success_count += 1,
//!         Err(e) => {
//!             eprintln!("Skipping molecule: {}", e);
//!             error_count += 1;
//!         }
//!     }
//! }
//! println!("Parsed {} molecules, {} errors", success_count, error_count);
//! ```
//!
//! ## Working with Properties
//!
//! SDF files can contain key-value properties in the data block. These are
//! stored as a `HashMap<String, String>` on the molecule.
//!
//! ### Getting Properties
//!
//! ```rust
//! use sdfrust::parse_sdf_string;
//!
//! let sdf = r#"aspirin
//!
//!
//!   1  0  0  0  0  0  0  0  0  0999 V2000
//!     0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
//! M  END
//! > <PUBCHEM_CID>
//! 2244
//!
//! > <MOLECULAR_WEIGHT>
//! 180.16
//!
//! $$$$
//! "#;
//!
//! let mol = parse_sdf_string(sdf).unwrap();
//!
//! // Get a single property
//! if let Some(cid) = mol.get_property("PUBCHEM_CID") {
//!     assert_eq!(cid, "2244");
//! }
//!
//! // Check if property exists
//! assert!(mol.properties.contains_key("MOLECULAR_WEIGHT"));
//! ```
//!
//! ### Setting Properties
//!
//! ```rust
//! use sdfrust::Molecule;
//!
//! let mut mol = Molecule::new("example");
//! mol.set_property("SMILES", "CCO");
//! mol.set_property("SOURCE", "generated");
//!
//! assert_eq!(mol.get_property("SMILES"), Some("CCO"));
//! ```
//!
//! ### Iterating Over Properties
//!
//! ```rust
//! use sdfrust::Molecule;
//!
//! let mut mol = Molecule::new("example");
//! mol.set_property("MW", "180.16");
//! mol.set_property("CID", "2244");
//!
//! for (key, value) in &mol.properties {
//!     println!("{}: {}", key, value);
//! }
//! ```
//!
//! ## Molecule Operations
//!
//! The `Molecule` struct provides many useful methods for working with
//! chemical structure data.
//!
//! ### Molecular Formula
//!
//! ```rust
//! use sdfrust::parse_sdf_string;
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
//! assert_eq!(mol.formula(), "H2O");
//! ```
//!
//! ### Geometric Center (Centroid)
//!
//! ```rust
//! use sdfrust::{Molecule, Atom};
//!
//! let mut mol = Molecule::new("example");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 2.0, 0.0, 0.0));
//!
//! let (cx, cy, cz) = mol.centroid().unwrap();
//! assert!((cx - 1.0).abs() < 1e-6);
//! assert!((cy - 0.0).abs() < 1e-6);
//! ```
//!
//! ### Bond Connectivity (Neighbors)
//!
//! ```rust
//! use sdfrust::{Molecule, Atom, Bond, BondOrder};
//!
//! let mut mol = Molecule::new("methane");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "H", 1.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "H", -1.0, 0.0, 0.0));
//! mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
//! mol.bonds.push(Bond::new(0, 2, BondOrder::Single));
//!
//! // Get neighbors of the carbon (index 0)
//! let neighbors = mol.neighbors(0);
//! assert_eq!(neighbors.len(), 2);
//! assert!(neighbors.contains(&1));
//! assert!(neighbors.contains(&2));
//! ```
//!
//! ### Element Counts
//!
//! ```rust
//! use sdfrust::{Molecule, Atom};
//!
//! let mut mol = Molecule::new("ethanol");
//! mol.atoms.push(Atom::new(0, "C", 0.0, 0.0, 0.0));
//! mol.atoms.push(Atom::new(1, "C", 1.5, 0.0, 0.0));
//! mol.atoms.push(Atom::new(2, "O", 2.5, 0.0, 0.0));
//! mol.atoms.push(Atom::new(3, "H", 0.0, 1.0, 0.0));
//!
//! let counts = mol.element_counts();
//! assert_eq!(counts.get("C"), Some(&2));
//! assert_eq!(counts.get("O"), Some(&1));
//! assert_eq!(counts.get("H"), Some(&1));
//! ```
//!
//! ### Other Useful Methods
//!
//! - `atom_count()` / `bond_count()` - Get counts
//! - `is_empty()` - Check if molecule has atoms
//! - `total_charge()` - Sum of formal charges
//! - `has_aromatic_bonds()` - Check for aromaticity
//! - `has_charges()` - Check for charged atoms
//! - `atoms_by_element("C")` - Filter atoms by element
//! - `bonds_by_order(BondOrder::Double)` - Filter bonds by type
//! - `translate(dx, dy, dz)` - Move molecule
//! - `center()` - Move centroid to origin
//!
//! ## Performance Tips
//!
//! ### Use Iterators for Large Files
//!
//! For files with thousands of molecules, use the iterator API to process
//! molecules one at a time without loading all into memory:
//!
//! ```rust,ignore
//! use sdfrust::iter_sdf_file;
//!
//! // Memory efficient - processes one molecule at a time
//! for result in iter_sdf_file("large_database.sdf")? {
//!     let mol = result?;
//!     // Process and discard
//! }
//!
//! // vs. loading all at once (uses more memory)
//! let all_molecules = parse_sdf_file_multi("large_database.sdf")?;
//! ```
//!
//! ### Release Builds for Benchmarks
//!
//! Parsing performance improves significantly with optimizations:
//!
//! ```bash
//! cargo build --release
//! cargo run --release --example benchmark
//! ```
//!
//! ### Streaming vs Load-All Tradeoffs
//!
//! | Approach | Memory | Speed | Use Case |
//! |----------|--------|-------|----------|
//! | `iter_sdf_file` | O(1) | Fast | Large files, filtering |
//! | `parse_sdf_file_multi` | O(n) | Fast | Need random access |
//! | `parse_sdf_string` | O(1) | Fastest | Single molecule |
//!
//! ## Format Notes
//!
//! ### Supported Formats
//!
//! - **SDF V2000**: Full support for reading and writing (up to 999 atoms/bonds)
//! - **MOL2 TRIPOS**: Full support for reading (MOLECULE, ATOM, BOND sections)
//!
//! ### SDF V3000
//!
//! SDF V3000 format is fully supported for both parsing and writing:
//!
//! ```rust
//! use sdfrust::{parse_sdf_auto_string, write_sdf_auto_string, SdfFormat};
//!
//! // V3000 content is automatically detected and parsed
//! let v3000_content = r#"test
//!
//!
//!   0  0  0     0  0            999 V3000
//! M  V30 BEGIN CTAB
//! M  V30 COUNTS 2 1 0 0 0
//! M  V30 BEGIN ATOM
//! M  V30 1 C 0.0000 0.0000 0.0000 0
//! M  V30 2 O 1.2000 0.0000 0.0000 0
//! M  V30 END ATOM
//! M  V30 BEGIN BOND
//! M  V30 1 2 1 2
//! M  V30 END BOND
//! M  V30 END CTAB
//! M  END
//! $$$$
//! "#;
//!
//! let mol = parse_sdf_auto_string(v3000_content).unwrap();
//! assert_eq!(mol.format_version, SdfFormat::V3000);
//! ```
//!
//! ### Format Detection
//!
//! The library uses file content to determine format:
//! - SDF V2000 files contain `V2000` in the counts line
//! - SDF V3000 files contain `V3000` in the counts line
//! - MOL2 files start with `@<TRIPOS>MOLECULE`

pub mod atom;
pub mod bond;
pub mod collection;
pub mod descriptors;
pub mod error;
pub mod featurize;
pub mod fingerprints;
#[cfg(feature = "geometry")]
pub mod geometry;
pub mod graph;
pub mod molecule;
pub mod parser;
pub mod sgroup;
pub mod stereogroup;
pub mod writer;

// Re-export main types
pub use atom::Atom;
pub use bond::{Bond, BondOrder, BondStereo};
pub use collection::{Collection, CollectionType};
pub use error::{Result, SdfError};
pub use molecule::{Molecule, SdfFormat};
pub use sgroup::{SGroup, SGroupType};
pub use stereogroup::{StereoGroup, StereoGroupType};

// Re-export parser functions
pub use parser::{
    SdfIterator, SdfParser, detect_sdf_format, iter_sdf_file, parse_sdf_auto_file,
    parse_sdf_auto_file_multi, parse_sdf_auto_string, parse_sdf_auto_string_multi, parse_sdf_file,
    parse_sdf_file_multi, parse_sdf_string, parse_sdf_string_multi,
};

// Re-export unified auto-detection functions and types
pub use parser::{
    AutoIterator, FileFormat, detect_format, iter_auto_file, parse_auto_file,
    parse_auto_file_multi, parse_auto_string, parse_auto_string_multi,
};

// Re-export V3000 parser functions
pub use parser::{
    SdfV3000Iterator, SdfV3000Parser, iter_sdf_v3000_file, parse_sdf_v3000_file,
    parse_sdf_v3000_file_multi, parse_sdf_v3000_string, parse_sdf_v3000_string_multi,
};

// Re-export MOL2 parser functions
pub use parser::{
    Mol2Iterator, Mol2Parser, iter_mol2_file, parse_mol2_file, parse_mol2_file_multi,
    parse_mol2_string, parse_mol2_string_multi,
};

// Re-export XYZ parser functions
pub use parser::{
    XyzIterator, XyzParser, iter_xyz_file, parse_xyz_file, parse_xyz_file_multi, parse_xyz_string,
    parse_xyz_string_multi,
};

// Re-export writer functions
pub use writer::{
    write_sdf, write_sdf_file, write_sdf_file_multi, write_sdf_multi, write_sdf_string,
};

// Re-export V3000 writer functions
pub use writer::{
    needs_v3000, write_sdf_auto, write_sdf_auto_file, write_sdf_auto_string, write_sdf_v3000,
    write_sdf_v3000_file, write_sdf_v3000_file_multi, write_sdf_v3000_multi,
    write_sdf_v3000_string,
};

// Re-export MOL2 writer functions
pub use writer::{
    write_mol2, write_mol2_file, write_mol2_file_multi, write_mol2_multi, write_mol2_string,
};

// Re-export bond inference functions
pub use descriptors::{BondInferenceConfig, infer_bonds, infer_bonds_with_config};

// Re-export graph module
pub use graph::AdjacencyList;

// Re-export ML featurization
pub use featurize::ogb::{OgbAtomFeatures, OgbBondFeatures, OgbGraphFeatures};

// Re-export fingerprints
pub use fingerprints::ecfp::EcfpFingerprint;

// Re-export descriptors module (access via sdfrust::descriptors::*)
// For direct access: use sdfrust::descriptors::{molecular_weight, exact_mass, ...}
