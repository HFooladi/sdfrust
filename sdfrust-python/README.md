# sdfrust - Python Bindings

Fast Rust-based SDF and MOL2 molecular structure file parser with Python bindings.

## Installation

### From source (requires Rust toolchain)

```bash
cd sdfrust-python
pip install maturin
maturin develop --features numpy
```

### Build wheel

```bash
maturin build --release --features numpy
pip install target/wheels/sdfrust-*.whl
```

## Quick Start

```python
import sdfrust

# Parse a single SDF file
mol = sdfrust.parse_sdf_file("molecule.sdf")
print(f"Name: {mol.name}")
print(f"Atoms: {mol.num_atoms}")
print(f"Formula: {mol.formula()}")
print(f"MW: {mol.molecular_weight():.2f}")

# Parse multiple molecules
mols = sdfrust.parse_sdf_file_multi("database.sdf")
for mol in mols:
    print(f"{mol.name}: {mol.num_atoms} atoms")

# Memory-efficient iteration over large files
for mol in sdfrust.iter_sdf_file("large_database.sdf"):
    print(f"{mol.name}: MW={mol.molecular_weight():.2f}")
```

## Supported Formats

- **SDF V2000**: Full support for reading and writing (up to 999 atoms/bonds)
- **SDF V3000**: Full support for reading and writing (unlimited atoms/bonds)
- **MOL2 TRIPOS**: Full support for reading

## API Reference

### Parsing Functions

#### SDF Files

```python
# Single molecule
mol = sdfrust.parse_sdf_file("file.sdf")      # V2000
mol = sdfrust.parse_sdf_auto_file("file.sdf") # Auto-detect V2000/V3000
mol = sdfrust.parse_sdf_v3000_file("file.sdf") # V3000 only

# Multiple molecules
mols = sdfrust.parse_sdf_file_multi("file.sdf")
mols = sdfrust.parse_sdf_auto_file_multi("file.sdf")

# From string
mol = sdfrust.parse_sdf_string(content)
mols = sdfrust.parse_sdf_string_multi(content)
```

#### MOL2 Files

```python
mol = sdfrust.parse_mol2_file("file.mol2")
mols = sdfrust.parse_mol2_file_multi("file.mol2")
mol = sdfrust.parse_mol2_string(content)
```

#### Iterators (Memory-Efficient)

```python
for mol in sdfrust.iter_sdf_file("large.sdf"):
    process(mol)

for mol in sdfrust.iter_mol2_file("large.mol2"):
    process(mol)
```

### Writing Functions

```python
# Single molecule
sdfrust.write_sdf_file(mol, "output.sdf")
sdfrust.write_sdf_auto_file(mol, "output.sdf")  # Auto V2000/V3000
sdf_string = sdfrust.write_sdf_string(mol)

# Multiple molecules
sdfrust.write_sdf_file_multi(mols, "output.sdf")
```

### Molecule Properties

```python
mol = sdfrust.parse_sdf_file("aspirin.sdf")

# Basic info
print(mol.name)           # Molecule name
print(mol.num_atoms)      # Number of atoms
print(mol.num_bonds)      # Number of bonds
print(mol.formula())      # Molecular formula

# Descriptors
print(mol.molecular_weight())    # Molecular weight
print(mol.exact_mass())          # Monoisotopic mass
print(mol.heavy_atom_count())    # Non-hydrogen atoms
print(mol.ring_count())          # Number of rings
print(mol.rotatable_bond_count()) # Rotatable bonds
print(mol.total_charge())        # Sum of formal charges

# Geometry
centroid = mol.centroid()        # (x, y, z) center
mol.translate(1.0, 0.0, 0.0)     # Move molecule
mol.center()                     # Center at origin

# Properties (from SDF data block)
cid = mol.get_property("PUBCHEM_CID")
mol.set_property("SOURCE", "generated")
```

### Atom Access

```python
# Iterate over atoms
for atom in mol.atoms:
    print(f"{atom.element} at ({atom.x}, {atom.y}, {atom.z})")

# Get specific atom
atom = mol.get_atom(0)
print(atom.element)
print(atom.formal_charge)
print(atom.coords())  # (x, y, z) tuple

# Filter atoms
carbons = mol.atoms_by_element("C")
neighbors = mol.neighbors(0)  # Atom indices bonded to atom 0
```

### Bond Access

```python
# Iterate over bonds
for bond in mol.bonds:
    print(f"{bond.atom1}-{bond.atom2}: {bond.order}")

# Filter bonds
double_bonds = mol.bonds_by_order(sdfrust.BondOrder.double())
aromatic = mol.has_aromatic_bonds()

# Bond properties
bond = mol.bonds[0]
print(bond.is_aromatic())
print(bond.contains_atom(0))
print(bond.other_atom(0))  # Other atom in bond
```

### NumPy Integration

```python
import numpy as np
import sdfrust

mol = sdfrust.parse_sdf_file("molecule.sdf")

# Get coordinates as NumPy array
coords = mol.get_coords_array()  # Shape: (N, 3)
print(coords.shape)

# Modify and set back
coords[:, 0] += 10.0  # Translate in x
mol.set_coords_array(coords)

# Get atomic numbers
atomic_nums = mol.get_atomic_numbers()  # Shape: (N,)
```

### Creating Molecules

```python
import sdfrust

# Create empty molecule
mol = sdfrust.Molecule("water")

# Add atoms
mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))

# Add bonds
mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
mol.add_bond(sdfrust.Bond(0, 2, sdfrust.BondOrder.single()))

# Write to file
sdfrust.write_sdf_file(mol, "water.sdf")
```

## Performance

sdfrust is implemented in Rust for maximum performance. Benchmarks show it is
significantly faster than pure Python parsers and comparable to C++ implementations.

For large files, use the iterator API (`iter_sdf_file`) to process molecules
one at a time without loading the entire file into memory.

## License

MIT License
