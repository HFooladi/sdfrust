# CLAUDE.md

Pure-Rust library for parsing/writing SDF (V2000/V3000), MOL2, and XYZ chemical structure files. Python bindings via PyO3 in `sdfrust-python/`.

## Build & Test

```bash
cargo build                    # Build
cargo test                     # Run all tests (260+ Rust + doc-tests)
cargo clippy                   # Lint
cargo build --workspace        # Include Python bindings crate

# Python bindings
cd sdfrust-python && source .venv/bin/activate
maturin develop --features numpy
pytest tests/ -v               # 41 Python tests

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
├── error.rs               # SdfError enum (17 variants)
├── atom.rs                # Atom struct (index, element, x/y/z, charge)
├── bond.rs                # Bond, BondOrder, BondStereo
├── molecule.rs            # Molecule container + SdfFormat enum
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
│   ├── elements.rs        # Periodic table data
│   ├── molecular.rs       # MW, exact mass, heavy atom count
│   └── topological.rs     # Ring count, rotatable bonds
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
- `Molecule`: `atom_count()`, `bond_count()`, `formula()`, `centroid()`, `neighbors(idx)`, `element_counts()`, `atoms_by_element(elem)`, `get_property(key)`, `set_property(key, val)`
- `Atom`: fields `index`, `element`, `x`, `y`, `z`, `formal_charge`
- `Bond`: fields `atom1`, `atom2`, `order` (BondOrder enum), `stereo`

## Development Workflows

**Adding a feature**: Create module in `src/` → add to `lib.rs` exports → write tests in `tests/` → update ROADMAP.md

**Adding a file format**: Create parser in `src/parser/<fmt>.rs` → map to `Molecule` → add integration tests with real files → add writer if needed → update ROADMAP.md

## Test Data

Test files in `tests/test_data/`: `aspirin.sdf`, `caffeine_pubchem.sdf`, `glucose.sdf`, `galactose.sdf`, `acetaminophen.sdf`, `methionine.sdf`, `methane.mol2`, `benzene.mol2`, `water.xyz`, `multi.xyz`

## Python Examples

Example scripts in `sdfrust-python/examples/` with PubChem drug data in `examples/data/`:

- `basic_usage.py` — Core API: parsing, writing, atoms, bonds, descriptors, NumPy
- `format_conversion.py` — Multi-format detection, XYZ parsing, SDF/MOL2 conversion, round-trips
- `batch_analysis.py` — Drug library processing: filtering, sorting, Lipinski analysis
- `geometry_analysis.py` — 3D geometry: distance matrices, RMSD, rotation, transforms (requires `--features geometry`)

## CI/CD

GitHub Actions: `.github/workflows/rust.yml`. Use `gh run list`, `gh run view <id>`, `gh pr list`.

## Format Gotchas

- **SDF V2000**: Fixed-width columns. Coordinates at positions 0-9, 10-19, 20-29. Element at 31-33. Bond atoms at 0-2, 3-5 (1-based!). Charge codes are inverted: 1=+3, 2=+2, 3=+1, 5=-1, 6=-2, 7=-3.
- **SDF V3000**: Variable-width, prefixed with `M  V30`. Supports >999 atoms, stereogroups, sgroups, extended bond types (hydrogen, ionic, coordination).
- **MOL2**: Section headers `@<TRIPOS>`. SYBYL atom types like `C.ar` — element extracted before `.`. Bond type `ar` for aromatic.
- **XYZ**: No bonds. Element can be symbol or atomic number. Case-normalized. Multi-molecule files are concatenated blocks.

## Roadmap

Phase 9.7 (XYZ Parser) complete. Next: Phase 10 (Shared Traits). See ROADMAP.md.
