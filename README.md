# sdfrust

[![Crates.io](https://img.shields.io/crates/v/sdfrust.svg)](https://crates.io/crates/sdfrust)
[![docs.rs](https://docs.rs/sdfrust/badge.svg)](https://docs.rs/sdfrust)
[![CI](https://github.com/hfooladi/sdfrust/actions/workflows/rust.yml/badge.svg)](https://github.com/hfooladi/sdfrust/actions/workflows/rust.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fast, pure-Rust parser for SDF (Structure Data File), MOL2, and XYZ chemical structure files, with Python bindings.

## Features

- **SDF V2000/V3000**: Full read/write support for both SDF format versions
- **MOL2 (TRIPOS)**: Full read/write support for MOL2 format
- **XYZ**: Read support for XYZ coordinate files (single and multi-molecule)
- **Gzip Support**: Transparent decompression of `.gz` files (optional `gzip` feature)
- **Python Bindings**: First-class Python API with NumPy integration
- **Streaming Parsing**: Memory-efficient iteration over large files
- **Molecular Descriptors**: MW, exact mass, ring count, rotatable bonds, and more
- **ML-Ready Featurization**: OGB-compatible GNN features, ECFP fingerprints, Gasteiger charges
- **Chemical Perception**: SSSR rings, aromaticity (Hückel 4n+2), hybridization, conjugation
- **3D GNN Support**: Cutoff-based neighbor lists, bond/dihedral angles for SchNet/DimeNet/GemNet
- **High Performance**: ~220,000 molecules/sec (4-7x faster than RDKit)
- **Real-World Validated**: 100% success rate on PDBbind 2024 dataset (27,670 ligand SDF files)

## Installation

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
sdfrust = "0.6"
```

### Python

```bash
# From source (requires Rust toolchain)
cd sdfrust-python
pip install maturin
maturin develop --features numpy
```

## Quick Start

### Rust

```rust
use sdfrust::{parse_sdf_file, parse_mol2_file, parse_xyz_file, parse_auto_file, write_sdf_file};

// Parse a single molecule (any format)
let mol = parse_sdf_file("molecule.sdf")?;
let mol = parse_mol2_file("molecule.mol2")?;
let mol = parse_xyz_file("coords.xyz")?;
let mol = parse_auto_file("unknown_format.sdf")?; // Auto-detect
println!("Name: {}", mol.name);
println!("Atoms: {}", mol.atom_count());
println!("Formula: {}", mol.formula());
println!("MW: {:.2}", mol.molecular_weight().unwrap());

// Parse multiple molecules
let molecules = parse_sdf_file_multi("database.sdf")?;

// Iterate over large files (memory efficient)
for result in iter_sdf_file("large_database.sdf")? {
    let mol = result?;
    // Process each molecule
}

// Write molecules
write_sdf_file("output.sdf", &mol)?;
```

### Python

```python
import sdfrust

# Parse molecules
mol = sdfrust.parse_sdf_file("molecule.sdf")
mol = sdfrust.parse_mol2_file("molecule.mol2")

# Access properties
print(f"Name: {mol.name}")
print(f"Atoms: {mol.num_atoms}")
print(f"Formula: {mol.formula()}")
print(f"MW: {mol.molecular_weight():.2f}")

# Molecular descriptors
print(f"Rings: {mol.ring_count()}")
print(f"Rotatable bonds: {mol.rotatable_bond_count()}")
print(f"Heavy atoms: {mol.heavy_atom_count()}")

# NumPy integration
import numpy as np
coords = mol.get_coords_array()  # (N, 3) array
atomic_nums = mol.get_atomic_numbers()  # (N,) array

# Iterate over large files
for mol in sdfrust.iter_sdf_file("large_database.sdf"):
    print(f"{mol.name}: MW={mol.molecular_weight():.2f}")

# Write molecules
sdfrust.write_sdf_file(mol, "output.sdf")
```

## Data Model

### Molecule

The main container for molecular data:

```rust
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub properties: HashMap<String, String>,
    pub format_version: SdfFormat,  // V2000 or V3000
}
```

### Atom

Represents an atom with 3D coordinates:

```rust
pub struct Atom {
    pub index: usize,
    pub element: String,
    pub x: f64, pub y: f64, pub z: f64,
    pub formal_charge: i8,
    pub mass_difference: i8,
    // ... additional fields
}
```

### Bond

Represents a bond between two atoms:

```rust
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
    pub stereo: BondStereo,
}
```

### BondOrder

```rust
pub enum BondOrder {
    Single, Double, Triple, Aromatic,
    SingleOrDouble, SingleOrAromatic, DoubleOrAromatic, Any,
    Coordination, Hydrogen,  // V3000 only
}
```

## Molecule Methods

```rust
// Atom/bond counts
mol.atom_count()
mol.bond_count()
mol.is_empty()

// Connectivity
mol.neighbors(atom_index)      // Get connected atom indices
mol.bonds_for_atom(atom_index) // Get bonds for an atom

// Properties
mol.get_property("MW")
mol.set_property("SMILES", "CCO")

// Chemistry
mol.formula()         // "C2H6O"
mol.total_charge()    // Sum of formal charges
mol.element_counts()  // HashMap of element counts
mol.has_aromatic_bonds()
mol.has_charges()

// Geometry
mol.centroid()        // Geometric center
mol.center()          // Move centroid to origin
mol.translate(x, y, z)

// Descriptors
mol.molecular_weight()     // IUPAC 2021 atomic weights
mol.exact_mass()           // Monoisotopic mass
mol.heavy_atom_count()     // Non-hydrogen atoms
mol.ring_count()           // Using Euler formula
mol.rotatable_bond_count() // RDKit-compatible definition
mol.is_atom_in_ring(idx)
mol.is_bond_in_ring(idx)

// Filtering
mol.atoms_by_element("C")
mol.bonds_by_order(BondOrder::Double)
```

## ML-Ready Features

sdfrust provides chemical perception and featurization for molecular ML — compute features in Rust, output NumPy arrays for PyTorch/JAX.

### OGB GNN Featurization

```rust
use sdfrust::featurize::ogb;

let graph = ogb::ogb_graph_features(&mol);
// graph.atom_features: [N, 9] - atomic_num, chirality, degree, charge, num_hs, radical, hybridization, aromatic, in_ring
// graph.bond_features: [2E, 3] - bond_type, stereo, conjugated
// graph.edge_src, graph.edge_dst: directed edge index
```

```python
# Python: direct NumPy arrays for PyTorch Geometric
atom_feats = mol.get_ogb_atom_features_array()  # np.array [N, 9]
bond_feats = mol.get_ogb_bond_features_array()  # np.array [E, 3]
graph = mol.ogb_graph_features()                 # dict with edge_src, edge_dst
```

### ECFP/Morgan Fingerprints

```rust
use sdfrust::fingerprints::ecfp;

let fp = ecfp::ecfp(&mol, 2, 2048);  // ECFP4, 2048 bits
let similarity = fp.tanimoto(&other_fp);
```

```python
fp = mol.ecfp(radius=2, n_bits=2048)       # List[bool]
fp_array = mol.get_ecfp_array()             # np.array (2048,)
sim = mol.tanimoto_similarity(other_mol)    # float
```

### Gasteiger Partial Charges

```rust
use sdfrust::descriptors::gasteiger_charges;

let charges = gasteiger_charges(&mol);  // Vec<f64>
```

```python
charges = mol.gasteiger_charges()           # List[float]
charges_array = mol.get_gasteiger_charges_array()  # np.array (N,)
```

### Chemical Perception

```rust
use sdfrust::descriptors::{sssr, all_aromatic_atoms, all_hybridizations, all_conjugated_bonds};

let rings = sssr(&mol);                         // SSSR ring perception
let aromatic = all_aromatic_atoms(&mol);         // Hückel 4n+2 aromaticity
let hybs = all_hybridizations(&mol);             // SP/SP2/SP3
let conjugated = all_conjugated_bonds(&mol);     // Bond conjugation
```

```python
rings = mol.sssr()                    # List of ring atom lists
aromatic = mol.all_aromatic_atoms()   # List[bool]
hybs = mol.all_hybridizations()       # List[str]
conj = mol.all_conjugated_bonds()     # List[bool]
```

### 3D GNN Features (geometry feature)

```rust
use sdfrust::geometry::{neighbor_list, all_bond_angles, all_dihedral_angles};

let nl = neighbor_list(&mol, 5.0);         // Cutoff-based neighbor list
let angles = all_bond_angles(&mol);         // Bond angles for DimeNet
let dihedrals = all_dihedral_angles(&mol);  // Dihedral angles for GemNet
```

```python
nl = mol.neighbor_list(cutoff=5.0)         # dict: edge_src, edge_dst, distances
angles = mol.all_bond_angles()              # dict: triplets, angles
dihedrals = mol.all_dihedral_angles()       # dict: quadruplets, angles
```

## Supported Formats

| Format | Read | Write |
|--------|------|-------|
| SDF V2000 | ✅ | ✅ |
| SDF V3000 | ✅ | ✅ |
| MOL2 (TRIPOS) | ✅ | ✅ |
| XYZ | ✅ | - |
| Gzip (`.gz`) | ✅ | - |

### Format Auto-Detection

```rust
// Automatically detect V2000 vs V3000
let mol = parse_sdf_auto_file("molecule.sdf")?;

// Automatically choose format based on molecule requirements
write_sdf_auto_file("output.sdf", &mol)?;
```

## Performance

The library is designed for high performance:

| Tool | Throughput | vs sdfrust |
|------|------------|------------|
| **sdfrust** | ~220,000 mol/s | baseline |
| RDKit | ~30,000-50,000 mol/s | 4-7x slower |
| Pure Python | ~3,000-5,000 mol/s | 40-70x slower |

- Streaming parser for memory-efficient processing of large files
- Minimal allocations during parsing
- Zero-copy where possible

## Python API Reference

### Parsing Functions

```python
# SDF (V2000/V3000)
mol = sdfrust.parse_sdf_file(path)
mol = sdfrust.parse_sdf_string(content)
mol = sdfrust.parse_sdf_auto_file(path)  # Auto-detect format
mols = sdfrust.parse_sdf_file_multi(path)

# MOL2
mol = sdfrust.parse_mol2_file(path)
mol = sdfrust.parse_mol2_string(content)
mols = sdfrust.parse_mol2_file_multi(path)

# Iterators (memory-efficient)
for mol in sdfrust.iter_sdf_file(path):
    process(mol)
```

### Writing Functions

```python
sdfrust.write_sdf_file(mol, path)
sdfrust.write_sdf_string(mol)
sdfrust.write_sdf_auto_file(mol, path)  # Auto-select V2000/V3000
sdfrust.write_sdf_file_multi(mols, path)
```

### ML Feature Methods (on Molecule)

```python
# Chemical perception
mol.sssr()                           # SSSR rings
mol.all_aromatic_atoms()             # Aromaticity
mol.all_hybridizations()             # Hybridization
mol.all_conjugated_bonds()           # Conjugation

# ML features
mol.ogb_atom_features()              # [N, 9] OGB features
mol.ogb_bond_features()              # [E, 3] OGB features
mol.ogb_graph_features()             # Full graph dict
mol.ecfp(radius=2, n_bits=2048)      # ECFP bit vector
mol.tanimoto_similarity(other)       # Tanimoto coefficient
mol.gasteiger_charges()              # Partial charges

# NumPy arrays (with numpy feature)
mol.get_ogb_atom_features_array()    # np.array [N, 9]
mol.get_ogb_bond_features_array()    # np.array [E, 3]
mol.get_ecfp_array()                 # np.array (n_bits,)
mol.get_gasteiger_charges_array()    # np.array (N,)

# 3D features (with geometry feature)
mol.neighbor_list(cutoff=5.0)        # Neighbor list dict
mol.bond_angle(i, j, k)             # Angle in radians
mol.dihedral_angle(i, j, k, l)      # Dihedral in radians
```

### Classes

- `Molecule`: Main molecular container
- `Atom`: Atom with coordinates
- `Bond`: Bond between atoms
- `BondOrder`: Bond type enum
- `BondStereo`: Stereochemistry enum
- `SdfFormat`: Format version enum

## Citation

If you use sdfrust in your research, please cite:

```bibtex
@software{fooladi2025sdfrust,
  author = {Fooladi, Hosein},
  title = {sdfrust: A fast, pure-Rust parser for SDF, MOL2, and XYZ chemical structure files},
  year = {2025},
  url = {https://github.com/hfooladi/sdfrust},
  license = {MIT}
}
```

## License

MIT License - Copyright (c) 2025-2026 Hosein Fooladi
