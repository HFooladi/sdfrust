# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Transparent Gzip Decompression** (`src/parser/compression.rs`)
  - Optional feature: enable with `--features gzip`
  - Automatic decompression of `.gz` files based on extension
  - `MaybeGzReader` enum for transparent handling of plain and gzipped files
  - Helper functions: `open_maybe_gz()`, `is_gzip_path()`, `read_maybe_gz_to_string()`
  - Works with all parsers: SDF V2000, V3000, MOL2, XYZ
  - Works with all file functions: `parse_*_file()`, `parse_*_file_multi()`, `iter_*_file()`
  - Case-insensitive extension matching (`.gz`, `.GZ`, `.Gz`)
  - `GzipNotEnabled` error with helpful message when feature is disabled

- **Python Gzip Bindings**
  - `gzip_enabled()` function to check if gzip support is compiled in
  - All file parsing functions transparently support `.gz` files
  - All iterators support gzipped files
  - 16 pytest tests for gzip functionality

- **XYZ Format Parser** (`src/parser/xyz.rs`)
  - Parse XYZ molecular coordinate files (coordinates only, no bonds)
  - `XyzParser<R>` and `XyzIterator<R>` for streaming parsing
  - Support for atomic numbers as element identifiers (1 → H, 6 → C, 8 → O)
  - Element symbol case normalization (ca → Ca)
  - Multi-molecule XYZ files (concatenated blocks)
  - `parse_xyz_string()`, `parse_xyz_file()`, `iter_xyz_file()` functions

- **Unified Auto-Detection for XYZ**
  - `FileFormat::Xyz` variant added to format detection
  - `detect_format()` now recognizes XYZ format
  - `parse_auto_string()`, `parse_auto_file()`, `iter_auto_file()` support XYZ

- **Python Bindings for XYZ**
  - `parse_xyz_file()`, `parse_xyz_string()` functions
  - `parse_xyz_file_multi()`, `parse_xyz_string_multi()` functions
  - `iter_xyz_file()` streaming iterator
  - `XyzIterator` class for memory-efficient iteration

- **Test Coverage**
  - 32 XYZ integration tests
  - 12 XYZ unit tests
  - 8 Python XYZ tests
  - Test data files: `water.xyz`, `multi.xyz`

## [0.4.0] - 2026-01-30

### Added

- **Python Bindings** (`sdfrust-python/` crate)
  - PyO3-based Python bindings with full API coverage
  - `PyAtom`, `PyBond`, `PyBondOrder`, `PyBondStereo` wrapper classes
  - `PyMolecule` with all properties and descriptor methods
  - SDF V2000/V3000 parsing functions (`parse_sdf_file`, `parse_sdf_string`, etc.)
  - SDF writing functions (`write_sdf_file`, `write_sdf_string`, etc.)
  - MOL2 parsing functions (`parse_mol2_file`, `parse_mol2_string`, etc.)
  - Memory-efficient streaming iterators (`iter_sdf_file`, `iter_mol2_file`)
  - NumPy integration for coordinate arrays (`get_coords_array`, `set_coords_array`)
  - Atomic number array support (`get_atomic_numbers`)
  - 33 pytest tests covering all functionality

- **Workspace Configuration**
  - Root `Cargo.toml` now defines workspace with `sdfrust` and `sdfrust-python` members
  - Shared dependencies and build configuration

### Python API Highlights

```python
import sdfrust

# Parse molecules
mol = sdfrust.parse_sdf_file("molecule.sdf")
mol = sdfrust.parse_mol2_file("molecule.mol2")

# Access properties
print(mol.name, mol.num_atoms, mol.formula())
print(mol.molecular_weight(), mol.ring_count())

# NumPy integration
coords = mol.get_coords_array()  # (N, 3) array

# Iterate over large files
for mol in sdfrust.iter_sdf_file("large.sdf"):
    process(mol)
```

## [0.3.0] - 2026-01-30

### Added

- **Basic Molecular Descriptors** (`src/descriptors/` module)
  - `molecular_weight()` - Calculate molecular weight using IUPAC 2021 atomic weights
  - `exact_mass()` - Calculate monoisotopic mass
  - `heavy_atom_count()` - Count non-hydrogen atoms
  - `bond_type_counts()` - Count bonds by type (single, double, aromatic, etc.)
  - `ring_count()` - Count rings using Euler characteristic formula
  - `ring_atoms()` / `ring_bonds()` - Identify atoms and bonds in rings
  - `rotatable_bond_count()` - Count rotatable bonds (RDKit-compatible definition)

- **Element Data Table** (`src/descriptors/elements.rs`)
  - ~30 common elements with atomic weights and monoisotopic masses
  - `get_element()`, `atomic_weight()`, `monoisotopic_mass()` functions
  - Support for deuterium (D) and tritium (T)

- **Convenience Methods on Molecule**
  - `mol.molecular_weight()`, `mol.exact_mass()`
  - `mol.heavy_atom_count()`, `mol.bond_type_counts()`
  - `mol.ring_count()`, `mol.rotatable_bond_count()`
  - `mol.is_atom_in_ring(idx)`, `mol.is_bond_in_ring(idx)`

- **Comprehensive Test Suite**
  - 54 descriptor integration tests
  - Validated against PubChem reference data (aspirin, caffeine, etc.)

## [0.2.0] - 2025-01-29

### Added

- **SDF V3000 Parser and Writer**
  - Full V3000 format parsing with atoms, bonds, stereogroups, SGroups, collections
  - V3000 format writing with proper `M  V30` prefix formatting
  - Automatic format detection (V2000 vs V3000)
  - Auto-format selection for writing based on molecule requirements
  - Extended bond types: Coordination (9) and Hydrogen (10)
  - New types: `StereoGroup`, `SGroup`, `Collection`
  - `detect_sdf_format()`, `parse_sdf_auto_string()`, `parse_sdf_auto_file()`
  - `write_sdf_auto()`, `needs_v3000()` methods

- **MOL2 Parser**
  - Full support for TRIPOS MOL2 format
  - `Mol2Parser<R>` and `Mol2Iterator<R>` for streaming parsing
  - Parse SYBYL atom types and MOL2 bond types (single, double, triple, aromatic)
  - `parse_mol2_string()`, `parse_mol2_file()`, `iter_mol2_file()` functions

- **Performance Benchmarking Suite**
  - Criterion benchmarks for SDF and MOL2 parsing
  - Python comparison benchmarks (vs RDKit and pure Python)
  - Results: ~220,000 mol/s (4-7x faster than RDKit, 40-70x faster than pure Python)

- **Comprehensive Testing**
  - 169 total tests (unit, integration, edge cases, round-trip)
  - Real-world test files from PubChem (aspirin, caffeine, glucose, etc.)
  - V3000-specific test data files

- **Documentation**
  - Expanded lib.rs documentation with examples and guides
  - CLAUDE.md for AI assistance
  - LICENSE (MIT) and CONTRIBUTORS.md

## [0.1.0] - 2025-01-25

### Added

- **SDF V2000 Parser**
  - Parse atom block (coordinates, element, charge, mass difference)
  - Parse bond block (atom indices, bond order, stereo)
  - Parse properties block (key-value pairs after M END)
  - Handle multi-molecule files ($$$$-delimited)
  - M CHG lines (formal charges) and M ISO lines (isotopes)
  - Memory-efficient iterator for large files

- **SDF V2000 Writer**
  - Write SDF V2000 format with property preservation
  - Write M CHG/M ISO lines for special atoms
  - Multi-molecule file output

- **Core Data Structures**
  - `Molecule`: Container with atoms, bonds, and properties
  - `Atom`: 3D coordinates, element, charge, mass difference
  - `Bond`: Atom indices, order, stereochemistry
  - `BondOrder` enum: Single, Double, Triple, Aromatic, etc.
  - `BondStereo` enum: None, Up, Down, Either

- **Molecule Methods**
  - `formula()`, `centroid()`, `translate()`, `center()`
  - `neighbors()`, `bonds_for_atom()`, `element_counts()`

[Unreleased]: https://github.com/HFooladi/sdfrust/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/HFooladi/sdfrust/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/HFooladi/sdfrust/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HFooladi/sdfrust/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HFooladi/sdfrust/releases/tag/v0.1.0
