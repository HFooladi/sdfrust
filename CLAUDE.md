# CLAUDE.md

Pure-Rust library for parsing/writing SDF (V2000/V3000), MOL2, and XYZ chemical structure files. Python bindings via PyO3 in `sdfrust-python/`.

## Build & Test

```bash
cargo build                                     # Build
cargo test                                       # Run all ~530 Rust tests (unit + integration + doc-tests)
cargo test --features geometry                   # Include geometry-gated tests
cargo clippy --workspace --features geometry     # Lint everything
cargo fmt --check                                # Check formatting
cargo build --workspace                          # Include Python bindings crate

# Python bindings
cd sdfrust-python && source .venv/bin/activate
maturin develop --features numpy
pytest tests/ -v                                 # ~346 Python tests

# PDBbind benchmark (~27k files, requires dataset)
PDBBIND_2024_DIR=/path/to/PDBbind_2024 cargo test --release pdbbind_benchmark -- --ignored --nocapture
```

## Code Conventions

- `thiserror` for error types; return `Result<T, SdfError>` from fallible operations
- Prefer iterators over index loops
- Use `BufRead` trait for parser input (not concrete types)
- All public items need doc comments
- Internal indices are 0-based (file formats use 1-based — convert at parse/write boundaries)
- Atoms/bonds own their data (no lifetimes)
- Properties stored as `HashMap<String, String>` on `Molecule`

## Architecture

```
src/
├── lib.rs                 # Public API re-exports
├── error.rs               # SdfError enum (18 variants)
├── atom.rs                # Atom struct (index, element, x/y/z, charge)
├── bond.rs                # Bond, BondOrder, BondStereo
├── molecule.rs            # Molecule container + SdfFormat enum
├── graph.rs               # AdjacencyList: pre-computed neighbor lookups, degree vectors
├── parser/
│   ├── sdf.rs             # V2000 parser + auto-detection + unified parse_auto_*
│   ├── sdf_v3000.rs       # V3000 parser
│   ├── mol2.rs            # MOL2 parser
│   └── xyz.rs             # XYZ parser
├── writer/
│   ├── sdf.rs             # V2000 writer
│   ├── sdf_v3000.rs       # V3000 writer + auto-format selection
│   └── mol2.rs            # MOL2 writer
├── descriptors/
│   ├── elements.rs        # Periodic table data + covalent radii
│   ├── molecular.rs       # MW, exact mass, heavy atom count
│   ├── topological.rs     # Ring count, rotatable bonds
│   ├── bond_inference.rs  # Infer bonds from 3D coordinates
│   ├── valence.rs         # Atom degree, implicit/total hydrogen count
│   ├── rings.rs           # SSSR ring perception (spanning-tree + GF(2))
│   ├── hybridization.rs   # SP/SP2/SP3/SP3D/SP3D2 from bond topology
│   ├── aromaticity.rs     # Hückel 4n+2 aromaticity perception
│   ├── conjugation.rs     # Conjugated bond detection
│   ├── chirality.rs       # CIP R/S chirality perception
│   └── gasteiger.rs       # Gasteiger-Marsili partial charges (PEOE)
├── geometry/              # Feature-gated: geometry = ["nalgebra"]
│   ├── neighbor_list.rs   # Cutoff-based neighbor list for 3D GNNs
│   └── angles.rs          # Bond angles + dihedral angles
├── featurize/
│   └── ogb.rs             # OGB-compatible GNN featurizer (9 atom + 3 bond features)
├── fingerprints/
│   └── ecfp.rs            # ECFP/Morgan fingerprints (Rogers & Hahn 2010)
├── sgroup.rs              # SGroup types (V3000)
├── stereogroup.rs         # Stereogroup types (V3000)
└── collection.rs          # Collection types (V3000)

sdfrust-python/src/        # PyO3 bindings (mirrors Rust API)
```

Each parser follows the same pattern: `Parser<R: BufRead>` for streaming + `Iterator<R>` for multi-molecule files. Top-level convenience functions: `parse_*_file()`, `parse_*_string()`, `parse_*_file_multi()`, `iter_*_file()`.

## Key API

- `parse_sdf_file(path)` / `parse_sdf_string(s)` — single molecule
- `parse_sdf_file_multi(path)` — all molecules into Vec
- `iter_sdf_file(path)` — streaming iterator (memory-efficient)
- `parse_auto_file(path)` — auto-detect format (SDF/MOL2/XYZ)
- `write_sdf_file(mol, path)` / `write_sdf_auto_file(mol, path)` — V2000 or auto V2000/V3000
- `infer_bonds(mol, tolerance)` / `infer_bonds_with_config(mol, config)` — infer single bonds from 3D coordinates
- `Molecule`: `atom_count()`, `bond_count()`, `formula()`, `centroid()`, `neighbors(idx)`, `element_counts()`, `atoms_by_element(elem)`, `get_property(key)`, `set_property(key, val)`, `infer_bonds(tolerance)`
- `Atom`: fields `index`, `element`, `x`, `y`, `z`, `formal_charge`
- `Bond`: fields `atom1`, `atom2`, `order` (BondOrder enum), `stereo`

### ML Feature API

- `AdjacencyList::from_molecule(mol)` — graph adjacency with O(1) lookups
- `descriptors::sssr(mol)` — SSSR ring perception
- `descriptors::atom_hybridization(mol, idx)` — SP/SP2/SP3 hybridization
- `descriptors::is_aromatic_atom(mol, idx)` / `all_aromatic_atoms(mol)` — aromaticity
- `descriptors::is_conjugated_bond(mol, idx)` / `all_conjugated_bonds(mol)` — conjugation
- `descriptors::gasteiger_charges(mol)` — partial charges (PEOE)
- `descriptors::atom_chirality(mol, idx)` / `all_chiralities(mol)` — CIP R/S chirality
- `featurize::ogb::ogb_atom_features(mol)` — [N, 9] OGB atom features (all 9 match RDKit)
- `featurize::ogb::ogb_bond_features(mol)` — [E, 3] OGB bond features
- `featurize::ogb::ogb_graph_features(mol)` — complete graph with directed edge index
- `fingerprints::ecfp::ecfp(mol, radius, n_bits)` — ECFP/Morgan fingerprint
- `geometry::neighbor_list(mol, cutoff)` — cutoff-based neighbor list (geometry feature)
- `geometry::bond_angle(mol, i, j, k)` / `dihedral_angle(mol, i, j, k, l)` — angles (geometry feature)

## Development Workflows

**Adding a feature**: Create module in `src/` → add to `lib.rs` exports → write tests in `tests/` → update ROADMAP.md

**Adding a file format**: Create parser in `src/parser/<fmt>.rs` → map to `Molecule` → add integration tests with real files → add writer if needed → update ROADMAP.md

## Test Data

Test files in `tests/test_data/`: `aspirin.sdf`, `caffeine_pubchem.sdf`, `glucose.sdf`, `galactose.sdf`, `acetaminophen.sdf`, `methionine.sdf`, `v3000_benzene.sdf`, `v3000_charged.sdf`, `v3000_methane.sdf`, `v3000_stereo.sdf`, `methane.mol2`, `benzene.mol2`, `water.xyz`, `multi.xyz`

Additional test molecules in `sdfrust-python/examples/data/`: `ibuprofen.sdf`, `dopamine.sdf`, `cholesterol.sdf`, `drug_library.sdf` (6 molecules)

## Python Examples

Example scripts in `sdfrust-python/examples/` with PubChem drug data in `examples/data/`:

- `basic_usage.py` — Core API: parsing, writing, atoms, bonds, descriptors, NumPy
- `format_conversion.py` — Multi-format detection, XYZ parsing, SDF/MOL2 conversion, round-trips
- `batch_analysis.py` — Drug library processing: filtering, sorting, Lipinski analysis
- `geometry_analysis.py` — 3D geometry: distance matrices, RMSD, rotation, transforms (requires `--features geometry`)
- `ml_features.py` — ML features: OGB featurization, ECFP fingerprints, Gasteiger charges, rings, aromaticity, NumPy arrays

## CI/CD

GitHub Actions: `.github/workflows/rust.yml`. Use `gh run list`, `gh run view <id>`, `gh pr list`.

## Gotchas

- **Python venv + Conda**: If `maturin develop` fails with linker errors, run `unset CONDA_PREFIX` before building
- **SDF wedge bonds**: Stereo field on single bonds encodes chirality (Up/Down), not E/Z — only encode bond stereo on double bonds

## Format Gotchas

- **SDF V2000**: Fixed-width columns. Coordinates at positions 0-9, 10-19, 20-29. Element at 31-33. Bond atoms at 0-2, 3-5 (1-based!). Charge codes are inverted: 1=+3, 2=+2, 3=+1, 5=-1, 6=-2, 7=-3.
- **SDF V3000**: Variable-width, prefixed with `M  V30`. Supports >999 atoms, stereogroups, sgroups, extended bond types (hydrogen, ionic, coordination).
- **MOL2**: Section headers `@<TRIPOS>`. SYBYL atom types like `C.ar` — element extracted before `.`. Bond type `ar` for aromatic.
- **XYZ**: No bonds. Element can be symbol or atomic number. Case-normalized. Multi-molecule files are concatenated blocks.

## Roadmap

v0.6.0 released: Phase 11 (ML Features, Tiers 0-4) complete including CIP chirality and bond inference. Next: Phase 11.12 (batch pipeline). See ROADMAP.md.
