# sdfrust Python Examples

This directory contains example scripts demonstrating the Python API.

## Setup

Before running examples, build and install the Python bindings:

```bash
cd sdfrust-python
uv venv .venv --python 3.11
source .venv/bin/activate
uv pip install maturin numpy pytest
maturin develop --features numpy
```

## Running Examples

```bash
# From sdfrust-python directory
python examples/basic_usage.py
```

## Examples Overview

### basic_usage.py

Comprehensive examples covering:

- **Parsing SDF files**: Load molecules from V2000 and V3000 SDF files
- **Parsing SDF strings**: Parse inline SDF content
- **Creating molecules**: Build molecules programmatically with atoms and bonds
- **Atom access**: Read atom properties, filter by element, calculate distances
- **Bond access**: Work with bond orders, stereochemistry, and connectivity
- **Molecular descriptors**: Calculate MW, exact mass, ring count, rotatable bonds
- **Geometry operations**: Centroid calculation, translation, centering
- **Properties**: Read/write SDF data block properties
- **Multi-molecule files**: Iterate over files with multiple molecules
- **Writing SDF**: Export molecules to V2000 or V3000 format
- **MOL2 parsing**: Read TRIPOS MOL2 format files
- **V3000 format**: Work with extended SDF format
- **NumPy integration**: Get/set coordinates as NumPy arrays
- **Ring detection**: Identify atoms and bonds in rings
- **Charged atoms**: Handle formal charges

## Quick Reference

```python
import sdfrust

# Parse files
mol = sdfrust.parse_sdf_file("molecule.sdf")
mol = sdfrust.parse_mol2_file("molecule.mol2")
mol = sdfrust.parse_sdf_auto_file("molecule.sdf")  # Auto-detect V2000/V3000

# Parse strings
mol = sdfrust.parse_sdf_string(sdf_content)
mol = sdfrust.parse_mol2_string(mol2_content)

# Create molecule
mol = sdfrust.Molecule("my_molecule")
mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))

# Basic properties
print(mol.name, mol.num_atoms, mol.num_bonds)
print(mol.formula())
print(mol.molecular_weight())

# Iterate large files efficiently
for mol in sdfrust.iter_sdf_file("large_library.sdf"):
    process(mol)

# Write files
sdfrust.write_sdf_file(mol, "output.sdf")
sdf_str = sdfrust.write_sdf_string(mol)

# NumPy integration
coords = mol.get_coords_array()  # shape (N, 3)
mol.set_coords_array(new_coords)
atomic_nums = mol.get_atomic_numbers()
```
