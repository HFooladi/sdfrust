# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.6.0] - 2026-02-19

### Added

- Graph adjacency infrastructure (`AdjacencyList`) with O(1) neighbor lookups and degree computation
- Atom valence and hydrogen count: `atom_degree()`, `implicit_hydrogen_count()`, `total_hydrogen_count()` with default valence table for ~20 elements
- SSSR ring perception via spanning-tree + GF(2) linear independence: `sssr()`, `ring_sizes()`, `smallest_ring_size()`
- Hybridization inference from bond topology: `atom_hybridization()`, `all_hybridizations()` (S/SP/SP2/SP3/SP3D/SP3D2)
- Aromaticity perception: two-stage (trust file annotations + Hückel 4n+2 rule on SSSR rings) with support for pyrrole, furan, thiophene
- Conjugation detection: aromatic bonds, SP2-SP2 single bonds, unsaturated bond chains
- OGB-compatible GNN featurizer: `ogb_atom_features()` [N, 9], `ogb_bond_features()` [E, 3], `ogb_graph_features()` with directed edge index
- ECFP/Morgan fingerprints: `ecfp()` bit vector, `ecfp_counts()` hash map, Tanimoto similarity — first pure-Rust implementation
- Gasteiger-Marsili partial charges: PEOE iterative algorithm with ~15 atom type parameters
- CIP chirality perception: `atom_chirality()`, `all_chiralities()` with R/S determination from 3D signed volume
- Bond inference from 3D coordinates: `infer_bonds()`, `infer_bonds_with_config()` using covalent radii
- Cutoff-based neighbor list for 3D GNNs: `neighbor_list()` with directed edges and distances (geometry feature)
- Bond angle and dihedral angle computation: `bond_angle()`, `dihedral_angle()`, `all_bond_angles()`, `all_dihedral_angles()` (geometry feature)
- Python bindings for all new ML features with NumPy array outputs (`get_ogb_atom_features_array()`, `get_ecfp_array()`, `get_gasteiger_charges_array()`)
- Python examples: format conversion, batch analysis, geometry operations, ML features
- Rust example: `examples/ml_features.rs`

## [0.5.0] - 2026-02-17

### Added

- Transparent gzip decompression for all file parsers (optional `gzip` feature)
- XYZ format parser with streaming support, atomic number input, and case normalization
- XYZ auto-detection in `parse_auto_*()` and `iter_auto_file()`
- Python bindings for XYZ parsing and gzip support
- PDBbind 2024 large-scale SDF benchmark test (27,670 files, 100% success rate, ~14k files/sec)

## [0.4.0] - 2026-01-30

### Added

- Python bindings via PyO3 (`sdfrust-python/` crate) with full API coverage
- NumPy integration for coordinate and atomic number arrays
- Streaming iterators for Python (`iter_sdf_file`, `iter_mol2_file`)
- Workspace configuration for multi-crate build

## [0.3.0] - 2026-01-30

### Added

- Molecular descriptors module: molecular weight, exact mass, heavy atom count, bond type counts
- Topological descriptors: ring count, ring atoms/bonds detection, rotatable bond count
- Element data table (~30 elements with IUPAC 2021 atomic weights)
- Convenience descriptor methods on `Molecule`

## [0.2.0] - 2025-01-29

### Added

- SDF V3000 parser and writer with stereogroups, SGroups, and collections
- Automatic V2000/V3000 format detection and auto-format writing
- MOL2 (TRIPOS) parser with streaming support
- Performance benchmarks (~220k mol/s, 4-7x faster than RDKit)

## [0.1.0] - 2025-01-25

### Added

- SDF V2000 parser with atom/bond/property blocks and M CHG/M ISO lines
- SDF V2000 writer with property preservation and multi-molecule output
- Core types: `Molecule`, `Atom`, `Bond`, `BondOrder`, `BondStereo`
- Molecule methods: `formula()`, `centroid()`, `translate()`, `neighbors()`, `element_counts()`
- Memory-efficient streaming iterator for large multi-molecule files

[Unreleased]: https://github.com/HFooladi/sdfrust/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/HFooladi/sdfrust/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/HFooladi/sdfrust/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/HFooladi/sdfrust/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/HFooladi/sdfrust/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HFooladi/sdfrust/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HFooladi/sdfrust/releases/tag/v0.1.0
