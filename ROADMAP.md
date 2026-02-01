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

**Status:** Implemented (integrated with Phase 6)

### Deliverables
- [x] Detect format from file content (V2000 vs V3000)
- [x] `detect_sdf_format()` function
- [x] `parse_sdf_auto_string()` / `parse_sdf_auto_file()` functions
- [x] `write_sdf_auto()` - auto-selects V2000/V3000 based on molecule needs
- [x] `needs_v3000()` - checks if molecule requires V3000 format
- [ ] Detect format from file extension (.sdf, .mol, .mol2) - deferred
- [ ] Gzip support (transparent decompression) - deferred

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
- [x] File I/O bindings for SDF and MOL2
- [x] NumPy array support for coordinates (`get_coords_array`, `set_coords_array`)
- [x] Atomic number array support (`get_atomic_numbers`)
- [x] Iterator support for large files (`iter_sdf_file`, `iter_mol2_file`)
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
    └── test_basic.py    # 33 pytest tests
```

### Test Coverage
- 33 pytest tests covering:
  - Version and module import
  - Atom creation and methods
  - Bond creation and methods
  - Molecule creation, atoms, bonds, properties
  - SDF string and file parsing
  - MOL2 string and file parsing
  - SDF writing
  - Iterators
  - Molecular descriptors
  - Geometry operations
  - NumPy coordinate arrays

### Python API
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

## Phase 11: Advanced Features

**Status:** Future

### Potential Features
- [ ] SMILES parsing/generation
- [ ] InChI generation
- [ ] Substructure search
- [ ] Fingerprint generation (ECFP, MACCS)
- [ ] 2D coordinate generation
- [ ] Stereochemistry perception
- [ ] Aromaticity perception
- [ ] Hydrogen addition/removal

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
| 1.0.0   | 10-11  | Stable API, advanced features |

---

## Contributing

Each phase should include:
1. Implementation
2. Unit tests
3. Integration tests
4. Documentation updates
5. CHANGELOG entry
