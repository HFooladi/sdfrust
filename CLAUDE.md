# CLAUDE.md - AI Development Guide for sdfrust

This document provides context for AI assistants working on this codebase.

## Project Overview

**sdfrust** is a pure-Rust library for parsing and writing SDF (Structure Data File) and MOL2 chemical structure files. It focuses on fast, safe file I/O for small molecule data.

## Architecture

### Core Modules

```
src/
├── lib.rs              # Public API exports
├── error.rs            # Error types (SdfError)
├── atom.rs             # Atom struct (coords, element, charge)
├── bond.rs             # Bond, BondOrder, BondStereo
├── molecule.rs         # Molecule container + methods
├── parser/
│   ├── mod.rs          # Parser exports
│   ├── sdf.rs          # SDF V2000 parser + iterator
│   ├── sdf_v3000.rs    # SDF V3000 parser + iterator
│   └── mol2.rs         # TRIPOS MOL2 parser + iterator
├── writer/
│   ├── mod.rs          # Writer exports
│   ├── sdf.rs          # SDF V2000 writer
│   └── sdf_v3000.rs    # SDF V3000 writer + auto-format
├── descriptors/
│   ├── mod.rs          # Descriptor exports
│   ├── elements.rs     # Element data table
│   ├── molecular.rs    # MW, exact mass, heavy atom count
│   └── topological.rs  # Ring count, rotatable bonds
├── sgroup.rs           # SGroup types (V3000)
├── stereogroup.rs      # Stereogroup types (V3000)
└── collection.rs       # Collection types (V3000)

sdfrust-python/         # Python bindings (PyO3)
├── src/
│   ├── lib.rs          # PyO3 module registration
│   ├── error.rs        # Error conversion
│   ├── atom.rs         # PyAtom wrapper
│   ├── bond.rs         # PyBond, PyBondOrder, PyBondStereo
│   ├── molecule.rs     # PyMolecule wrapper + NumPy support
│   ├── parsing.rs      # Parsing functions
│   ├── writing.rs      # Writing functions
│   └── iterators.rs    # Iterator wrappers
├── python/sdfrust/     # Python package
└── tests/              # pytest tests
```

### Key Types

- `Molecule`: Main container with atoms, bonds, and properties
- `Atom`: Represents an atom with 3D coordinates, element, charge
- `Bond`: Connects two atoms with order (single/double/triple/aromatic) and stereo
- `BondOrder`: Enum for bond types (1-10 as per SDF spec, including V3000 types)
- `SdfFormat`: Enum for V2000/V3000 format
- `SdfParser<R>`: Streaming SDF V2000 parser from any `BufRead`
- `SdfV3000Parser<R>`: Streaming SDF V3000 parser
- `SdfIterator<R>`: Iterator over molecules in multi-molecule SDF files
- `Mol2Parser<R>`: Streaming MOL2 parser from any `BufRead`
- `Mol2Iterator<R>`: Iterator over molecules in multi-molecule MOL2 files
- `StereoGroup`, `SGroup`, `Collection`: V3000-specific structures

### Design Patterns

1. **Streaming Parsing**: Uses `BufRead` trait for memory-efficient large file handling
2. **Zero-based Indexing**: Internal indices are 0-based (SDF uses 1-based)
3. **Properties as HashMap**: SDF data block stored as key-value pairs
4. **Owned Strings**: Atoms/bonds own their data (no lifetimes)

## Build Commands

```bash
# Rust library
cargo build              # Build library
cargo test               # Run all tests
cargo test pubchem       # Run PubChem validation tests
cargo doc --open         # Generate documentation
cargo clippy             # Run linter
cargo build --workspace  # Build all crates including Python bindings

# Python bindings
cd sdfrust-python
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install maturin numpy pytest
maturin develop --features numpy  # Build and install in dev mode
pytest tests/ -v                   # Run Python tests
```

## Test Structure

```
tests/                         # Rust integration tests
├── sdf_parsing_tests.rs       # Basic SDF V2000 functionality
├── pubchem_tests.rs           # Real PubChem molecules
├── edge_case_tests.rs         # Edge case handling
├── mol2_tests.rs              # MOL2 integration tests
├── v3000_tests.rs             # SDF V3000 tests
├── descriptor_tests.rs        # Molecular descriptor tests
└── test_data/
    ├── aspirin.sdf            # PubChem CID 2244
    ├── caffeine_pubchem.sdf   # PubChem CID 2519
    ├── glucose.sdf            # PubChem CID 5988 (sucrose)
    ├── galactose.sdf          # PubChem CID 5793
    ├── acetaminophen.sdf      # PubChem CID 1983
    ├── methionine.sdf         # PubChem CID 6137
    ├── methane.mol2           # Simple MOL2 test file
    └── benzene.mol2           # Aromatic MOL2 test file

sdfrust-python/tests/          # Python tests
└── test_basic.py              # 33 pytest tests
```

Total: 200+ Rust tests (unit + integration) + 33 Python tests

## SDF Format Reference

### V2000 Structure
```
Molecule Name           <- Line 1: name (any string)
  Program  Timestamp    <- Line 2: usually blank or program info
Comment line            <- Line 3: comment
aaabbb...V2000          <- Line 4: counts (aaa=atoms, bbb=bonds)
x y z elem ...          <- Atom block (one line per atom)
111222ttt...            <- Bond block (111=atom1, 222=atom2, ttt=type)
M  CHG  n ...           <- Property block
M  END                  <- End of mol
> <PROPERTY_NAME>       <- Data block
value
                        <- Blank line between properties
$$$$                    <- End of record
```

### Key Parsing Details

1. **Coordinates**: Positions 0-9, 10-19, 20-29 (10 chars each, 4 decimals)
2. **Element**: Positions 31-33 (3 chars, left-aligned)
3. **Charge Code**: Position 36-38 (0=none, 1=+3, 2=+2, 3=+1, 5=-1, 6=-2, 7=-3)
4. **Bond Atoms**: Positions 0-2, 3-5 (1-based indices!)
5. **Bond Type**: Position 6-8 (1=single, 2=double, 3=triple, 4=aromatic)

## Common Development Tasks

### Adding a New Feature

1. Create feature module in `src/`
2. Add to `lib.rs` exports
3. Write unit tests in `tests/`
4. Update ROADMAP.md

### Adding New File Format

1. Create parser in `src/parser/<format>.rs`
2. Map to `Molecule` struct
3. Add tests with real-world files
4. Update ROADMAP.md

### Running Specific Tests

```bash
cargo test test_aspirin       # Tests containing "aspirin"
cargo test pubchem --nocapture # PubChem tests with output
cargo test -- --test-threads=1 # Sequential execution
```

## GitHub CI/CD

GitHub Actions workflows are configured in `.github/workflows/rust.yml`. Use `gh` to interact with GitHub:

```bash
gh run list                    # List recent workflow runs
gh run view <run_id>           # View details of a specific run
gh run view <run_id> --log     # View full logs
gh run watch <run_id>          # Watch a run until completion
gh pr list                     # List pull requests
gh pr view <pr_number>         # View PR details
```

## Code Style

- Use `thiserror` for error types
- Return `Result<T, SdfError>` from fallible operations
- Prefer iterators over index loops
- Use `BufRead` trait for input (not concrete types)
- All public items need documentation

## Related Projects

- **pdbrust**: Sister project for PDB/mmCIF files (same author)
- **chemcore**: Another Rust cheminformatics library (proof-of-concept)
- **rdkit-rs**: RDKit bindings (requires C++)

## MOL2 Format Reference

### Structure
```
@<TRIPOS>MOLECULE         <- Section header
molecule_name             <- Molecule name
num_atoms num_bonds ...   <- Counts line
mol_type                  <- SMALL, BIOPOLYMER, etc.
charge_type               <- NO_CHARGES, USER_CHARGES, etc.

@<TRIPOS>ATOM             <- Atom section
atom_id name x y z type [subst_id subst_name charge]

@<TRIPOS>BOND             <- Bond section
bond_id atom1 atom2 type  <- type: 1, 2, 3, ar, am
```

### Key Parsing Details

1. **Section Headers**: Start with `@<TRIPOS>`
2. **Atom Types**: SYBYL types like "C.ar", "N.pl3" - element extracted before "."
3. **Bond Types**: "1" (single), "2" (double), "3" (triple), "ar" (aromatic)
4. **Atom Indices**: 1-based in file, converted to 0-based internally

## Roadmap

See [ROADMAP.md](ROADMAP.md) for detailed development phases.

Current status: Phase 9 (Python Bindings) complete. Next: Phase 10 (Shared Traits).
