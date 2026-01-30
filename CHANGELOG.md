# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/HFooladi/sdfrust/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/HFooladi/sdfrust/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HFooladi/sdfrust/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HFooladi/sdfrust/releases/tag/v0.1.0
