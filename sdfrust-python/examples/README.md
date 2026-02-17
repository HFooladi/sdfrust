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
python examples/format_conversion.py
python examples/batch_analysis.py

# Geometry examples require the geometry feature
maturin develop --features numpy,geometry
python examples/geometry_analysis.py
```

## Examples Overview

### basic_usage.py

Core functionality covering:

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

### format_conversion.py

Multi-format detection, parsing, and conversion:

- **Format detection**: Identify SDF V2000, V3000, and MOL2 from content (`detect_format()`)
- **Auto-detection parsing**: Parse any format with one function (`parse_auto_file()`)
- **XYZ parsing**: Single and multi-molecule XYZ files, streaming with `iter_xyz_file()`
- **SDF → MOL2 conversion**: Convert molecules with round-trip verification
- **MOL2 → SDF conversion**: Convert molecules with round-trip verification
- **Batch conversion**: Write multi-molecule SDF and MOL2 files (`write_sdf_file_multi()`, `write_mol2_file_multi()`)
- **Gzip support**: Check runtime feature availability (`gzip_enabled()`)

### batch_analysis.py

Compound library processing with real drug molecules:

- **Library loading**: Parse multi-molecule SDF with summary table
- **Streaming iteration**: Memory-efficient processing with `iter_sdf_file()`
- **MW filtering**: Filter by molecular weight range, write subsets
- **Descriptor sorting**: Rank molecules by MW and heavy atom count
- **Element composition**: Compare C/H/O/N/S counts across molecules
- **Bond profiles**: Compare single/double/aromatic bond distributions
- **Lipinski analysis**: Apply Rule of Five criteria using SDF properties
- **Finding extremes**: Identify molecules with most atoms, highest MW, most rings

### geometry_analysis.py

3D geometry operations (requires `--features geometry`):

- **Distance matrix**: Pairwise distances, closest/farthest atom pairs
- **Bond length analysis**: Extract actual bond lengths by bond type
- **RMSD comparison**: Compare identical, translated, and perturbed molecules
- **Rotation**: Axis-angle rotation with distance preservation verification
- **Combined transforms**: Rotation matrix + translation
- **NumPy analysis**: Bounding box, spatial extent, radius of gyration
- **Conformer comparison**: Perturb coordinates, center, compare RMSD

### Example Data

The `data/` directory contains real drug molecules downloaded from PubChem:

- `ibuprofen.sdf` — NSAID (CID 3672, 33 atoms)
- `dopamine.sdf` — Neurotransmitter (CID 681, 22 atoms)
- `cholesterol.sdf` — Steroid (CID 5997, 74 atoms)
- `drug_library.sdf` — Combined multi-molecule file (6 drugs)

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
