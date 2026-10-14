# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.8.0] - 2026-10-14

Chemistry correctness release. Every perceived property was compared against
RDKit 2025.09 and the `ogb` package on 4,553 molecules (NCI first 5k, EGFR and
CDK2 ligands, heavy-atom-only and explicit-H; 113k atoms) with
`compare_tools.py --rdkit-data`. Agreement before → after:

| Feature | 0.7.0 | 0.8.0 |
|---|---|---|
| Molecular weight / exact mass / formula | 54% | 100% |
| Aromatic atoms / bonds | 97.8% / 97.2% | 100% / 100% |
| Hybridization | 98.8% | >99.99% |
| Conjugated bonds | 98.4% | >99.99% |
| Gasteiger charges | 46% | 100% |
| Rotatable bonds | 68% | 100% |
| SSSR ring sizes | 99.8% | 100% |
| OGB bond stereo | 99.4% | 100% |
| ECFP bits identical to RDKit Morgan | ~0% | 100% |

### Changed (breaking: values change for the same input)

- `molecular_weight`, `exact_mass` and `Molecule::formula` count implicit
  hydrogens, so heavy-atom-only files give the full molecule. They also use
  isotope masses for labelled atoms, and exact mass corrects for the electron
  mass of ions. Formulas use RDKit's style: C, H, then alphabetical, with a
  charge suffix (`C2H3O2-`). A molecule built from atoms without bonds now
  gets implicit H on every atom; add the bonds, or set `Atom::no_implicit`.
- Implicit hydrogens follow RDKit's valence model: a charged atom takes the
  valences of its isoelectronic neutral element (N⁺ like C, B⁻² like N).
  Radicals and the molfile valence field are honoured.
- Protonated/quaternary N⁺ and sulfoxide/sulfonium S can be tetrahedral
  stereocenters (the lone pair counts as the fourth substituent).
- Aromaticity uses per-atom π-electron counts and fused-ring unions (RDKit's
  model). Kekulé naphthalene, indole, azulene, purines, pyridones and
  tropylium are recognized; sulfoxide rings are not.
- Hybridization comes from σ bonds plus lone pairs, and conjugation follows
  RDKit's rules. `atom_hybridization` and `all_hybridizations` now agree.
- Gasteiger charges are a port of RDKit's implementation: same parameter
  table, 12 iterations, implicit H as virtual atoms, H χ⁺ = 20.02, and formal
  charge split over conjugated atoms. Atoms without parameters keep their
  formal charge.
- `rotatable_bond_count` uses RDKit's default (strict) definition, which
  excludes amide/ester C–X bonds, CF3/CCl3/CBr3, tert-butyl and bonds next
  to triple bonds. Aspirin now has 2 rotatable bonds, not 3.
- `sssr` returns a true minimum cycle basis (Horton) with deterministic order.
- ECFP/Morgan fingerprints are bit-identical to RDKit's
  `GetMorganGenerator(radius, fpSize)` with default options.
- Chiral tags (CW/CCW) are relative to bond order, as in RDKit. Their docs no
  longer call them R/S.
- OGB bond stereo is perceived from coordinates: Z = 1, E = 2, unknown = 5.
- The element table covers all 118 elements, with isotope masses. MOL2 and
  XYZ atoms are flagged `no_implicit`.

Also checked on 2,000 PDBbind 2024 ligands (3D, explicit H; 123k atoms):
every property agrees with RDKit at ≥99.9%. All 27,670 PDBbind 2024 ligand
files still parse.
- The only ECFP or formal-charge differences are two N-oxides that RDKit
  rewrites while sanitizing (N=O becomes [n+][O-]); sdfrust keeps the charges
  as written in the file.
- 0.08% of chiral tags differ, on RDKit's ring-stereo centers (dependent
  stereocenters across rings), which aren't perceived yet.

### Added

- `descriptors::stereo::double_bond_stereo`: E/Z from coordinates.
- `gasteiger_charges_with_hydrogens`: also returns implicit-H charges.
- `element_by_atomic_number`, `isotope_mass`; `ElementData::outer_electrons`
  and `ElementData::valences`.
- `Atom::no_implicit`, also exposed in Python.
- `scripts/gen_element_table.py`, which regenerates the periodic table from RDKit.
- `compare_tools.py --rdkit-data/--smiles/--extra` for larger corpora, and
  exact ECFP bit comparison.

### Packaging

- Python wheels include the `geometry` feature (RMSD, transforms, neighbor
  lists, angles), so the documented geometry API works from PyPI.
- macOS and Windows wheels are built for Python 3.9–3.13; 0.6.0 shipped only
  3.12 wheels on those platforms.
- The Python sdist installs from source again. 0.6.0's failed because the
  crate's `[[bench]]` targets were excluded from the package. CI now builds
  and imports the sdist.

### Fixed

- The V2000 parser and writer handle stereo code 3 ("cis/trans unknown") on
  double bonds.

## [0.7.0] - not published

These changes were never released separately; they ship in 0.8.0.

Correctness release. Every change below was found by a differential comparison
against RDKit 2025.09, OpenBabel 3.1 and the `ogb` package
(`sdfrust-python/benchmarks/compare_tools.py`); each has a regression test.

### Changed (breaking)

- **OGB features now use the real `ogb.utils.features` encoding.** Values were
  previously raw values on a different scale, which is incompatible with OGB's
  `AtomEncoder`/`BondEncoder`. Models trained on 0.6 features must be retrained.
  - atomic number is `Z - 1` (misc = 118)
  - hybridization is `[SP, SP2, SP3, SP3D, SP3D2, misc]` (S → misc)
  - degree is total degree including implicit H
  - radicals are counted as electrons
  - bond stereo uses `[NONE, Z, E, CIS, TRANS, ANY]` (wedge flags are no longer read as E/Z)
  - out-of-range values go to the `misc` bucket
- `Hybridization::to_ogb_index()` returns OGB indices (SP = 0 … SP3D2 = 4, S/Other = 5).
- MOL2: bond type `am` (amide) is now `Single`, not `Aromatic`.
- MOL2: partial charges are stored in the new `Atom::partial_charge` instead of
  being rounded into `formal_charge`. `N.4` atoms get formal charge +1.
- Isotopes from `M  ISO` / V3000 `MASS=` are stored as absolute mass numbers in the
  new `Atom::isotope: Option<u16>`. Previously they went into `mass_difference: i8`,
  which also corrupted masses ≥ 128. `mass_difference` now only holds the V2000
  atom-block column.
- `write_sdf` (V2000) returns `SdfError::Unwritable` for molecules it cannot
  represent: more than 999 atoms/bonds, coordinates beyond ±9999.9999, or long atom
  symbols. It previously wrote a corrupt file. `write_sdf_auto` falls back to V3000
  in these cases.
- Writers return `SdfError::InvalidAtomIndex` for bonds or groups that reference
  missing atoms. They previously panicked or wrote invalid files.
- `MaybeGzReader::Gzip` wraps `MultiGzDecoder`.

### Fixed

- `parse_sdf_*` / `iter_sdf_file` silently returned an empty molecule for V3000
  records. V3000 records are now detected per record, so mixed V2000/V3000 files
  parse correctly.
- A trailing blank line after the final `$$$$` made `*_multi` fail on the whole file.
- Iterators (SDF V2000/V3000, MOL2) now resume at the next record after a malformed
  one, instead of stopping, as the crate docs described.
- Panics on non-ASCII text in V2000 fixed-width lines. Invalid UTF-8 (e.g. Latin-1
  property values) is now replaced instead of aborting the file.
- Huge declared atom counts no longer abort with an allocation failure.
- A missing `M  END` no longer swallows the following record(s).
- V2000 charges: all nonzero charges are written to `M  CHG`. The parser lets
  `M  CHG`/`M  RAD` supersede atom-block charges, as the CTfile spec requires.
  Previously a molecule with a +1 and a +4 charge lost the +1 when read by RDKit.
- Radicals: atom-block charge code 4 and `M  RAD` are parsed and `M  RAD` is
  written. Radicals no longer force V3000 output.
- `M  CHG`/`M  ISO` were written using `atom.index` instead of the atom's position.
- Data items now end at a blank line. Property output is sorted by name, so it is
  deterministic.
- MOL2: non-contiguous atom ids are resolved through an id map; element case is
  normalized (`CL` → `Cl`); `NO_CHARGES` is honored.
- XYZ: atomic numbers 1–118 are supported (previously ~25 elements).
- `ecfp(n_bits = 0)` no longer panics (returns an empty fingerprint; Python raises
  `ValueError`). Python `Molecule.add_bond` raises `IndexError` for bad atom indices.
- Gzip: multi-member (concatenated / bgzip) files are read completely.

### Added

- `Atom::isotope`, `Atom::partial_charge`, `Atom::radical_electrons()`. Python
  `Atom` exposes `isotope`, `partial_charge`, `radical` and `radical_electrons`.
- `SdfParser::skip_to_record_end()`, `SdfV3000Parser::skip_to_record_end()`.
- `sdfrust-python/benchmarks/compare_tools.py`: a differential comparison harness
  against RDKit, OpenBabel and OGB, with Markdown/JSON reports.
- `tests/regression_tests.rs` and `sdfrust-python/tests/test_regressions.py`.
- The Python OGB validation tests compare against `ogb.utils.features` itself, and
  CI installs RDKit and OGB so they run.

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

[Unreleased]: https://github.com/HFooladi/sdfrust/compare/v0.8.0...HEAD
[0.8.0]: https://github.com/HFooladi/sdfrust/compare/v0.6.0...v0.8.0
[0.6.0]: https://github.com/HFooladi/sdfrust/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/HFooladi/sdfrust/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/HFooladi/sdfrust/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/HFooladi/sdfrust/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HFooladi/sdfrust/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HFooladi/sdfrust/releases/tag/v0.1.0
