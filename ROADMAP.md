# sdfrust Roadmap

This document outlines the development phases for sdfrust, a pure-Rust library for parsing chemical structure files.

## Phase 1: SDF V2000 Parser (Core) ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Parse SDF V2000 atom block (coordinates, element, charge, mass diff)
- [x] Parse SDF V2000 bond block (atom indices, bond order, stereo)
- [x] Parse properties block (key-value pairs after M END)
- [x] Handle multi-molecule files ($$$$-delimited)
- [x] M CHG lines (formal charges)
- [x] M ISO lines (isotopes)
- [x] Basic validation and error handling
- [x] Memory-efficient iterator for large files

### Data Model
- [x] `Molecule` struct with atoms, bonds, properties
- [x] `Atom` struct with coordinates, element, charge
- [x] `Bond` struct with order and stereochemistry
- [x] `BondOrder` enum (Single, Double, Triple, Aromatic, etc.)
- [x] `BondStereo` enum (None, Up, Down, Either)

---

## Phase 2: SDF V2000 Writer ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Write SDF V2000 format
- [x] Preserve properties on round-trip
- [x] Write M CHG/M ISO lines for special atoms
- [x] Multi-molecule file output

---

## Phase 3: Comprehensive Testing & Validation ✅ COMPLETE

**Status:** Complete

### Deliverables
- [x] Unit tests for parsing/writing
- [x] Real-world file tests (6 PubChem molecules)
- [x] Edge case tests (41 tests for empty molecules, special chars, etc.)
- [x] Round-trip validation suite
- [x] Performance benchmarks (~220K molecules/sec)
- [ ] Comparison with RDKit/OpenBabel output (deferred)

### Test Coverage Achieved
- [x] 129 test cases (target was 100+)
- [x] Real-world files from PubChem
- [x] Edge cases: empty molecules, large coordinates, charges, stereo, multi-molecule

### Performance (release build)
- Parse rate: ~220,000 molecules/sec
- Round-trip: ~16,000 molecules/sec
- Property lookup: ~55 ns/call
- Centroid calculation: ~263 ns/call

---

## Phase 4: MOL2 Parser ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Parse TRIPOS MOL2 format
  - [x] @<TRIPOS>MOLECULE record
  - [x] @<TRIPOS>ATOM record
  - [x] @<TRIPOS>BOND record
  - [x] @<TRIPOS>SUBSTRUCTURE record (skipped, optional)
  - [x] @<TRIPOS>COMMENT record (skipped, optional)
- [x] Map to same `Molecule` structure
- [x] Partial charges converted to formal charges (rounded)
- [x] Multi-molecule MOL2 files
- [x] Memory-efficient iterator for large files

### Test Coverage
- [x] 6 unit tests in parser module
- [x] 17 integration tests with real MOL2 files
- [x] Bond types: single, double, triple, aromatic
- [x] Edge cases: charged molecules, extra sections, atom types

### Data Model Notes
- SYBYL atom types (e.g., "C.ar", "N.pl3") parsed to extract element
- Partial charges rounded to formal charges
- Future: add `atom_type: Option<String>` and `partial_charge: Option<f64>` if needed

---

## Phase 5: Performance Benchmarking & Comparison ✅ COMPLETE

**Status:** Implemented and tested

**Motivation:** One of the primary reasons for implementing in Rust is performance. This must be validated with rigorous comparisons against established tools.

### Deliverables
- [x] Create comprehensive benchmark suite (`benches/` directory)
- [x] Criterion benchmarks for SDF parsing
- [x] Criterion benchmarks for MOL2 parsing
- [x] Criterion benchmarks for round-trip operations
- [x] Benchmark against Python/RDKit
  - [x] Parse varying molecule counts (10, 100, 1K, 10K)
  - [x] Measure memory usage with psutil
  - [x] Compare throughput (molecules/second)
- [x] Benchmark against pure Python implementations
  - [x] Simple line-by-line SDF parser
  - [x] Show Rust advantage clearly
- [x] Document results in BENCHMARK_RESULTS.md
- [ ] Automated benchmark CI (track regressions) - deferred to CI setup
- [ ] Benchmark against OpenBabel (C++) - deferred

### Benchmark Suite Structure
```
benches/
├── sdf_parse_benchmark.rs    # Criterion SDF benchmarks
├── mol2_parse_benchmark.rs   # Criterion MOL2 benchmarks
├── roundtrip_benchmark.rs    # Parse + write cycle benchmarks
├── comparison/
│   ├── benchmark_rdkit.py    # RDKit comparison
│   ├── benchmark_pure_python.py # Pure Python baseline
│   ├── generate_report.py    # Generate markdown report
│   ├── run_all.sh            # Master benchmark script
│   └── requirements.txt      # Python dependencies
└── data/
    ├── README.md             # Dataset acquisition instructions
    └── generate_synthetic.py # Synthetic SDF generator
```

### Benchmark Groups

| File | Groups | Purpose |
|------|--------|---------|
| `sdf_parse_benchmark.rs` | `parse_single`, `parse_multi`, `parse_iterator`, `parse_real_files` | SDF parsing performance |
| `mol2_parse_benchmark.rs` | `mol2_parse_single`, `mol2_parse_multi`, `mol2_parse_iterator` | MOL2 parsing performance |
| `roundtrip_benchmark.rs` | `roundtrip_single`, `write_only`, `molecule_operations`, `property_operations` | Full cycle benchmarks |

### Results Summary

| Tool | Throughput | vs sdfrust |
|------|------------|------------|
| **sdfrust** | ~220,000 mol/s | baseline |
| RDKit | ~30,000-50,000 mol/s | 4-7x slower |
| Pure Python | ~3,000-5,000 mol/s | 40-70x slower |

### Usage

```bash
# Run Criterion benchmarks
cargo bench

# Run full comparison suite
cd benches/comparison && ./run_all.sh

# View HTML reports
open target/criterion/report/index.html
```

---

## Phase 6: SDF V3000 Parser & Writer ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Parse V3000 counts line (M V30 BEGIN CTAB)
- [x] Parse V3000 atom block (M V30 BEGIN ATOM)
- [x] Parse V3000 bond block (M V30 BEGIN BOND)
- [x] Handle extended features:
  - [x] Extended bond types (Coordination=9, Hydrogen=10)
  - [x] R-group labels
  - [x] Enhanced stereochemistry (StereoGroups)
  - [x] Atom-to-atom mapping
  - [x] Radical state
  - [x] SGroups (superatoms, polymers)
  - [x] Collections (atom lists, R-groups)
- [x] V3000 writer
- [x] Automatic format detection (V2000 vs V3000)
- [x] Auto-format selection for writing (based on molecule needs)

### New Types
- `SdfFormat` enum: V2000, V3000
- `StereoGroup`, `StereoGroupType`: Enhanced stereochemistry
- `SGroup`, `SGroupType`: Superatoms, polymers, etc.
- `Collection`, `CollectionType`: Atom lists, R-groups

### Extended Atom/Bond Fields
- Atoms: v3000_id, atom_atom_mapping, rgroup_label, radical
- Bonds: v3000_id, reacting_center, Coordination, Hydrogen bond types

### Test Coverage
- 31 V3000-specific tests
- Parsing: basic molecules, charges, radicals, aromatic bonds, stereo
- Writing: simple molecules, charges, properties, round-trip
- Edge cases: empty molecules, atoms-only, multi-molecule files
- Format detection and auto-selection

---

## Phase 7: Format Auto-Detection ✅ COMPLETE

**Status:** Implemented (integrated with Phase 6, extended in Phase 9.7)

### Deliverables
- [x] Detect format from file content (V2000 vs V3000 vs MOL2 vs XYZ)
- [x] `detect_sdf_format()` function (SDF V2000/V3000 only)
- [x] `detect_format()` function (all formats: SDF V2000, V3000, MOL2, XYZ)
- [x] `parse_sdf_auto_string()` / `parse_sdf_auto_file()` functions
- [x] `parse_auto_string()` / `parse_auto_file()` functions (all formats)
- [x] `write_sdf_auto()` - auto-selects V2000/V3000 based on molecule needs
- [x] `needs_v3000()` - checks if molecule requires V3000 format
- [ ] Detect format from file extension (.sdf, .mol, .mol2, .xyz) - deferred
- [x] Gzip support (transparent decompression) - see Phase 9.8

---

## Phase 8: Basic Descriptors ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Molecular weight calculation (IUPAC 2021 atomic weights)
- [x] Exact mass calculation (monoisotopic masses)
- [x] Atom count by element (via `element_counts()`)
- [x] Bond count by type (`bond_type_counts()`)
- [x] Heavy atom count
- [x] Rotatable bond count (RDKit-compatible SMARTS definition)
- [x] Ring detection (DFS-based cycle detection)

### Module Structure
```
src/descriptors/
├── mod.rs           # Module exports
├── elements.rs      # Element data table (~30 common elements)
├── molecular.rs     # molecular_weight, exact_mass, heavy_atom_count
└── topological.rs   # ring_count, ring_atoms, ring_bonds, rotatable_bond_count
```

### Test Coverage
- 54 integration tests in `tests/descriptor_tests.rs`
- 39 unit tests in descriptor modules
- Validated against PubChem reference data (aspirin, caffeine, glucose, etc.)

### API
Descriptor functions are available via:
- `sdfrust::descriptors::*` module functions
- Convenience methods on `Molecule` struct (e.g., `mol.molecular_weight()`)

---

## Phase 9: Python Bindings ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] PyO3 module setup with Maturin
- [x] `PyMolecule` wrapper class with all properties and methods
- [x] `PyAtom`, `PyBond`, `PyBondOrder`, `PyBondStereo` wrappers
- [x] `PySdfFormat` for V2000/V3000 handling
- [x] File I/O bindings for SDF, MOL2, and XYZ
- [x] NumPy array support for coordinates (`get_coords_array`, `set_coords_array`)
- [x] Atomic number array support (`get_atomic_numbers`)
- [x] Iterator support for large files (`iter_sdf_file`, `iter_mol2_file`, `iter_xyz_file`)
- [x] All molecular descriptors exposed (MW, ring count, etc.)
- [x] Maturin build configuration with workspace integration
- [ ] PyPI package publication (pending)

### Module Structure
```
sdfrust-python/
├── Cargo.toml           # PyO3 + numpy dependencies
├── pyproject.toml       # Maturin configuration
├── src/
│   ├── lib.rs           # Module registration
│   ├── error.rs         # SdfError → Python exception mapping
│   ├── atom.rs          # PyAtom wrapper
│   ├── bond.rs          # PyBond, PyBondOrder, PyBondStereo
│   ├── molecule.rs      # PyMolecule + NumPy support
│   ├── parsing.rs       # Parsing functions
│   ├── writing.rs       # Writing functions
│   └── iterators.rs     # Iterator wrappers
├── python/sdfrust/      # Python package
│   ├── __init__.py      # Re-exports
│   └── py.typed         # PEP 561 marker
└── tests/
    └── test_basic.py    # 41 pytest tests
```

### Test Coverage
- 41 pytest tests covering:
  - Version and module import
  - Atom creation and methods
  - Bond creation and methods
  - Molecule creation, atoms, bonds, properties
  - SDF string and file parsing
  - MOL2 string and file parsing
  - XYZ string and file parsing
  - SDF writing
  - Iterators (SDF, MOL2, XYZ)
  - Molecular descriptors
  - Geometry operations
  - NumPy coordinate arrays

### Python API
```python
import sdfrust

# Parse molecules (multiple formats)
mol = sdfrust.parse_sdf_file("molecule.sdf")
mol = sdfrust.parse_mol2_file("molecule.mol2")
mol = sdfrust.parse_xyz_file("molecule.xyz")
mol = sdfrust.parse_auto_file("any_format.xyz")  # Auto-detect

# Access properties
print(mol.name, mol.num_atoms, mol.formula())
print(mol.molecular_weight(), mol.ring_count())

# NumPy integration
coords = mol.get_coords_array()  # (N, 3) array

# Iterate over large files
for mol in sdfrust.iter_sdf_file("large.sdf"):
    process(mol)
```

---

## Phase 9.5: Geometry Module

**Status:** Complete

### Deliverables
- [x] Feature-gated geometry module (`geometry = ["nalgebra"]`)
- [x] Distance matrix calculation
- [x] RMSD calculation (without alignment)
- [x] 3D rotation operations (rotate by axis/angle)
- [x] Apply rotation matrix and transformation
- [x] Python bindings for geometry operations
- [x] Integration tests for geometry functions

### Module Structure
```
src/geometry/
├── mod.rs           # Module exports and Molecule extension methods
├── transform.rs     # Coordinate transformations (rotation, apply_transform)
├── distance.rs      # Distance calculations (distance_matrix)
└── rmsd.rs          # RMSD calculation
```

### API
```rust
// Rust API (with geometry feature)
use sdfrust::Molecule;

let mut mol = parse_sdf_file("molecule.sdf")?;

// Distance matrix
let matrix = mol.distance_matrix();

// Rotation (90° around Z-axis)
mol.rotate([0.0, 0.0, 1.0], std::f64::consts::PI / 2.0);

// RMSD
let other = parse_sdf_file("other.sdf")?;
let rmsd = mol.rmsd_to(&other)?;
```

---

## Phase 9.7: XYZ Format Parser ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Parse XYZ molecular coordinate format
  - [x] Atom count line (first line)
  - [x] Comment/title line (second line, used as molecule name)
  - [x] Atom lines (element x y z, whitespace-separated)
- [x] Multi-molecule XYZ files (concatenated blocks)
- [x] Memory-efficient iterator for large files
- [x] Automatic format detection (added to `detect_format()`)
- [x] Full Python bindings

### Special Features
- [x] Atomic numbers as element identifiers (1 → H, 6 → C, 8 → O, etc.)
- [x] Element symbol case normalization (ca → Ca, CL → Cl)
- [x] Extra columns after x/y/z are ignored
- [x] Blank lines between molecules handled gracefully
- [x] Scientific notation coordinates supported

### API
```rust
use sdfrust::{parse_xyz_file, parse_xyz_string, iter_xyz_file};

// Parse single molecule
let mol = parse_xyz_file("water.xyz")?;
let mol = parse_xyz_string(xyz_content)?;

// Parse multiple molecules
let mols = parse_xyz_file_multi("trajectory.xyz")?;

// Stream large files
for result in iter_xyz_file("large.xyz")? {
    let mol = result?;
    process(mol);
}

// Auto-detection works with XYZ
let mol = parse_auto_file("molecule.xyz")?;  // Detects XYZ format
```

### Python API
```python
import sdfrust

# Parse XYZ files
mol = sdfrust.parse_xyz_file("water.xyz")
mol = sdfrust.parse_xyz_string(xyz_content)

# Multi-molecule
mols = sdfrust.parse_xyz_file_multi("trajectory.xyz")

# Iterate over large files
for mol in sdfrust.iter_xyz_file("large.xyz"):
    process(mol)

# Auto-detection
mol = sdfrust.parse_auto_file("molecule.xyz")
fmt = sdfrust.detect_format(content)  # Returns "xyz"
```

### Test Coverage
- 12 unit tests in parser module
- 32 integration tests in `tests/xyz_tests.rs`
- 8 Python tests in `test_basic.py`
- Edge cases: atomic numbers, case normalization, scientific notation, blank lines

### Notes
- XYZ format contains no bond information; `mol.bonds` will be empty
- Molecule name is taken from the comment line (line 2)
- Format detection checks: first line is integer, third line has element + 3 floats

---

## Phase 9.8: Transparent Gzip Decompression ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Optional gzip feature (`--features gzip`)
- [x] Automatic decompression of `.gz` files based on extension
- [x] `MaybeGzReader` enum for transparent handling (no dynamic dispatch)
- [x] Helper functions: `open_maybe_gz()`, `is_gzip_path()`, `read_maybe_gz_to_string()`
- [x] Works with all parsers: SDF V2000, V3000, MOL2, XYZ
- [x] Works with all file functions: `parse_*_file()`, `parse_*_file_multi()`, `iter_*_file()`
- [x] Case-insensitive extension matching (`.gz`, `.GZ`, `.Gz`)
- [x] `GzipNotEnabled` error with helpful message when feature is disabled
- [x] Full Python bindings with `gzip_enabled()` function

### Module Structure
```
src/parser/
└── compression.rs     # MaybeGzReader enum and helper functions
```

### API
```rust
// Rust API (with gzip feature)
use sdfrust::{parse_sdf_file, iter_sdf_file};

// Transparently handles gzipped files
let mol = parse_sdf_file("molecule.sdf.gz")?;

// Iterators also work with gzipped files
for result in iter_sdf_file("large.sdf.gz")? {
    let mol = result?;
    process(mol);
}

// Auto-detection works with gzipped files
let mol = parse_auto_file("molecule.mol2.gz")?;

// Check gzip path detection
use sdfrust::parser::is_gzip_path;
assert!(is_gzip_path("test.sdf.gz"));
```

### Python API
```python
import sdfrust

# Check if gzip support is available
if sdfrust.gzip_enabled():
    # All file functions transparently support .gz files
    mol = sdfrust.parse_sdf_file("molecule.sdf.gz")

    # Iterators work too
    for mol in sdfrust.iter_sdf_file("large.sdf.gz"):
        process(mol)
```

### Test Coverage
- 24 Rust integration tests in `tests/gzip_tests.rs`
- 16 Python tests in `tests/test_gzip.py`
- Edge cases: empty files, case variations, plain files with gzip feature

### Build Options
```bash
# Build without gzip (default, smaller binary)
cargo build

# Build with gzip support
cargo build --features gzip

# Python bindings with gzip
maturin develop --features gzip
```

### Design Notes
- Uses `flate2` crate for gzip decompression
- `MaybeGzReader` enum avoids `Box<dyn BufRead>` overhead
- Extension-based detection (not magic bytes) for simplicity
- Read-only: writing gzip files not supported (use external tools)

---

## Phase 9.9: Bond Inference from 3D Coordinates ✅ COMPLETE

**Status:** Implemented and tested

### Deliverables
- [x] Covalent radii data table (Cordero et al. 2008) with `covalent_radius()` lookup
- [x] Bond inference algorithm: distance ≤ sum of covalent radii + tolerance
- [x] `infer_bonds()` simple API with optional tolerance
- [x] `infer_bonds_with_config()` full config API (tolerance, clear_existing_bonds)
- [x] `BondInferenceConfig` struct for fine-grained control
- [x] `Molecule::infer_bonds()` convenience method
- [x] Python bindings (`mol.infer_bonds()`)
- [x] `BondInferenceError` variant for unknown elements
- [ ] Phase B: Bond order assignment from valence constraints (future)

### Module Structure
```
src/descriptors/
└── bond_inference.rs    # Core algorithm, config, Molecule impl
src/descriptors/
└── elements.rs          # Added covalent_radius() + static table
```

### API
```rust
use sdfrust::{parse_xyz_file, infer_bonds, infer_bonds_with_config, BondInferenceConfig};

// Simple API
let mut mol = parse_xyz_file("water.xyz")?;
infer_bonds(&mut mol, None)?;           // default tolerance (0.45 A)
infer_bonds(&mut mol, Some(0.3))?;      // custom tolerance

// Convenience method
mol.infer_bonds(None)?;

// Full config
let config = BondInferenceConfig {
    tolerance: 0.3,
    clear_existing_bonds: false,
    ..Default::default()
};
infer_bonds_with_config(&mut mol, &config)?;
```

### Python API
```python
import sdfrust

mol = sdfrust.parse_xyz_file("water.xyz")
mol.infer_bonds()             # default tolerance
mol.infer_bonds(tolerance=0.3)  # custom tolerance
print(mol.num_bonds)          # 2
```

### Test Coverage
- 8 unit tests in `bond_inference.rs`
- 15 integration tests in `tests/bond_inference_tests.rs`
- Molecules tested: water, methane, H2, ethanol, CO2, benzene
- Edge cases: empty, single atom, unknown element, overlapping atoms, distant atoms
- Tolerance effects, existing bond clearing

### Notes
- All inferred bonds are single bonds (Phase A only)
- Phase B (bond order assignment using valence constraints) deferred to Phase 11
- Default tolerance of 0.45 A matches xyz2mol and Open Babel
- Minimum distance threshold (0.01 A) prevents bonding overlapping atoms

---

## Phase 10: Shared Traits (mol-core)

**Status:** Planned

### Deliverables
- [ ] Create `mol-core` crate with trait definitions
- [ ] `MolecularStructure` trait
- [ ] `AtomLike` trait
- [ ] `BondLike` trait
- [ ] Update sdfrust to implement traits
- [ ] Update pdbrust to implement traits
- [ ] Enable cross-format operations

---

## Phase 11: ML-Ready Chemical Perception & Featurization ✅ COMPLETE (Tiers 0-3)

**Status:** Tiers 0-3 complete; Tier 4 (CIP chirality, batch pipeline) planned

**Motivation:** Make sdfrust the data preprocessing layer for molecular ML — handling fast I/O, feature computation, and tensor output so PyTorch/JAX models can consume data directly, bypassing SMILES-based pipelines.

### Phase 11.0 — Graph Adjacency Infrastructure ✅

**Module:** `src/graph.rs`

- [x] `AdjacencyList` struct: pre-computed neighbor lookups, degree vectors, O(1) access
- [x] `from_molecule()`, `neighbors()`, `neighbor_atoms()`, `degree()`, `heavy_degree()`
- [x] Foundation for all subsequent phases

### Phase 11.1 — Atom Degree & Implicit Hydrogen Count ✅

**Module:** `src/descriptors/valence.rs`

- [x] Default valence table: maps (element, charge) → typical valence (~20 entries)
- [x] `atom_degree()`, `bond_order_sum()`
- [x] `implicit_hydrogen_count()`, `total_hydrogen_count()`
- [x] `all_atom_degrees()`, `all_implicit_hydrogen_counts()`, `all_total_hydrogen_counts()`
- [x] **Unlocks 2 of 9 OGB atom features** (degree, num_hs)

### Phase 11.2 — Neighbor List with Cutoff ✅

**Module:** `src/geometry/neighbor_list.rs` (under `geometry` feature)

- [x] `neighbor_list(mol, cutoff) → NeighborList` with edge_src, edge_dst, distances
- [x] `neighbor_list_with_self_loops()` variant
- [x] Directed pairs (both i→j and j→i) matching PyTorch Geometric format
- [x] **Enables SchNet/DimeNet/PaiNN/GemNet workflows**

### Phase 11.3 — SSSR Ring Perception ✅

**Module:** `src/descriptors/rings.rs`

- [x] Smallest Set of Smallest Rings via spanning-tree + BFS + GF(2) independence
- [x] `sssr(mol) → Vec<Ring>` with actual ring atom and bond indices
- [x] `ring_sizes()`, `smallest_ring_size()`, `is_in_ring_of_size()`
- [x] Prerequisite for aromaticity perception

### Phase 11.4 — Bond Angles & Dihedral Angles ✅

**Module:** `src/geometry/angles.rs` (under `geometry` feature)

- [x] `bond_angle(mol, i, j, k) → f64` (radians)
- [x] `dihedral_angle(mol, i, j, k, l) → f64` (radians)
- [x] `all_bond_angles()` / `all_dihedral_angles()` — enumerated from bonded paths
- [x] Returns `(triplet_indices, angles)` / `(quadruplet_indices, angles)` for DimeNet/GemNet

### Phase 11.5 — Hybridization ✅

**Module:** `src/descriptors/hybridization.rs`

- [x] `enum Hybridization { S, SP, SP2, SP3, SP3D, SP3D2, Other }`
- [x] Inferred from bond order sum + neighbor count
- [x] `to_ogb_index()` for direct OGB feature encoding
- [x] **Unlocks 1 OGB feature** + prerequisite for conjugation/Gasteiger

### Phase 11.6 — Aromaticity Perception ✅

**Module:** `src/descriptors/aromaticity.rs`

- [x] Two-stage: (1) trust `BondOrder::Aromatic` from file, (2) Hückel 4n+2 on SSSR rings
- [x] Pi-electron counting for C, N, O, S, Se, P
- [x] `is_aromatic_atom()`, `is_aromatic_bond()`, `all_aromatic_atoms()`, `all_aromatic_bonds()`
- [x] **Unlocks 1 OGB feature** (is_aromatic)

### Phase 11.7 — Conjugation Detection ✅

**Module:** `src/descriptors/conjugation.rs`

- [x] Aromatic bonds → conjugated
- [x] Single bond between two SP2 atoms → conjugated
- [x] Double/triple bond adjacent to unsaturation → conjugated
- [x] **Completes all 3 OGB bond features** (bond_type, stereo, is_conjugated)

### Phase 11.8 — OGB-Compatible GNN Featurizer ✅

**Module:** `src/featurize/ogb.rs`

- [x] `ogb_atom_features(mol) → [N, 9]` integer matrix
- [x] `ogb_bond_features(mol) → [E, 3]` integer matrix
- [x] `ogb_graph_features(mol) → OgbGraphFeatures` with directed edge index
- [x] Matches OGB `AtomEncoder`/`BondEncoder` conventions exactly
- [x] **Replaces entire RDKit preprocessing pipeline**

### Phase 11.9 — ECFP/Morgan Fingerprints ✅

**Module:** `src/fingerprints/ecfp.rs`

- [x] Rogers & Hahn (2010) algorithm: initial atom invariant → neighbor collection → fold to bits
- [x] `ecfp(mol, radius, n_bits) → EcfpFingerprint` with bit vector
- [x] `ecfp_counts(mol, radius) → EcfpCountFingerprint` with hash → count map
- [x] Tanimoto similarity, density, on_bits
- [x] **First pure-Rust ECFP implementation**

### Phase 11.10 — Gasteiger Partial Charges ✅

**Module:** `src/descriptors/gasteiger.rs`

- [x] PEOE algorithm: iterative charge equalization
- [x] Electronegativity parameters for ~15 atom types × hybridization states
- [x] `gasteiger_charges(mol)` (6 iterations, 0.5 damping)
- [x] `gasteiger_charges_with_params(mol, max_iter, damping)`

### Phase 11 Python Bindings ✅

- [x] All features exposed on `Molecule` class
- [x] NumPy array outputs: `get_ogb_atom_features_array()`, `get_ogb_bond_features_array()`, `get_gasteiger_charges_array()`, `get_ecfp_array()`
- [x] OGB graph features as dict matching PyTorch Geometric format
- [x] ECFP fingerprints with Tanimoto similarity
- [x] Neighbor list and angle computation (geometry feature)

### Phase 11 Test Coverage

- 585+ tests passing (unit + integration + doc-tests)
- Benzene, pyrrole, furan, cyclopentane, naphthalene, cubane tested for SSSR
- Kekulized and aromatic form aromaticity detection
- Butadiene conjugation chain
- Water H-O-H angle validation
- Zero clippy warnings

### Module Structure
```
src/
├── graph.rs                    # AdjacencyList, degree (Phase 11.0)
├── descriptors/
│   ├── valence.rs              # Default valence, implicit H (Phase 11.1)
│   ├── rings.rs                # SSSR ring perception (Phase 11.3)
│   ├── hybridization.rs        # SP/SP2/SP3 (Phase 11.5)
│   ├── aromaticity.rs          # Hückel-rule perception (Phase 11.6)
│   ├── conjugation.rs          # Conjugated bond detection (Phase 11.7)
│   └── gasteiger.rs            # Partial charges (Phase 11.10)
├── geometry/
│   ├── neighbor_list.rs        # Cutoff-based neighbor list (Phase 11.2)
│   └── angles.rs               # Bond + dihedral angles (Phase 11.4)
├── featurize/
│   ├── mod.rs
│   └── ogb.rs                  # OGB feature vectors (Phase 11.8)
└── fingerprints/
    ├── mod.rs
    └── ecfp.rs                 # Morgan/ECFP (Phase 11.9)
```

### Remaining Phases (Future)

- [ ] Phase 11.11: CIP chirality (R/S labels from CIP priority rules)
- [ ] Phase 11.12: Parallel batch pipeline (optional `rayon` dependency)
- [ ] Bond order assignment from valence constraints (Phase B of bond inference)
- [ ] SMILES parsing/generation
- [ ] Substructure search

---

## Quality Standards

Following pdbrust conventions:

### Code Quality
- Comprehensive error handling with `thiserror`
- Zero unsafe code
- Clippy clean (`-D warnings`)
- Rustfmt formatted
- Documentation for all public items

### Testing
- Unit tests for all modules
- Integration tests with real files
- Property-based tests where applicable
- Performance benchmarks
- CI/CD with GitHub Actions

### Documentation
- README with examples
- Inline documentation
- CHANGELOG maintenance
- CLAUDE.md for AI assistance

---

## Version Milestones

| Version | Phases | Description |
|---------|--------|-------------|
| 0.1.0   | 1-2    | SDF V2000 read/write ✅ |
| 0.2.0   | 3-7    | Testing, MOL2, benchmarks, SDF V3000 ✅ |
| 0.3.0   | 8      | Basic descriptors ✅ |
| 0.4.0   | 9      | Python bindings ✅ |
| 0.5.0   | 9.7-9.9 | XYZ parser, gzip support, bond inference ✅ |
| 0.6.0   | 11.0-11.10 | ML features: graph adjacency, valence, SSSR, aromaticity, hybridization, conjugation, OGB featurizer, ECFP, Gasteiger, neighbor list, angles ✅ |
| 1.0.0   | 10-11.12 | Stable API, CIP chirality, batch pipeline |

---

## Contributing

Each phase should include:
1. Implementation
2. Unit tests
3. Integration tests
4. Documentation updates
5. CHANGELOG entry
