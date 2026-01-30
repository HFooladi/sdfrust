#!/usr/bin/env python3
"""
Basic Usage Examples for sdfrust Python Bindings

This script demonstrates the core functionality of the sdfrust library
for parsing, analyzing, and writing molecular structure files.

Run this script from the sdfrust-python directory after building:
    cd sdfrust-python
    maturin develop --features numpy
    python examples/basic_usage.py
"""

import os
import tempfile

import sdfrust

# Path to test data
TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")


def example_version():
    """Check library version."""
    print("=" * 60)
    print("Library Version")
    print("=" * 60)
    print(f"sdfrust version: {sdfrust.__version__}")
    print()


def example_parse_sdf_file():
    """Parse a single molecule from an SDF file."""
    print("=" * 60)
    print("Parsing SDF File")
    print("=" * 60)

    # Parse aspirin from PubChem
    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    print(f"Molecule name: {mol.name}")
    print(f"Number of atoms: {mol.num_atoms}")
    print(f"Number of bonds: {mol.num_bonds}")
    print(f"Molecular formula: {mol.formula()}")
    print(f"Format version: {mol.format_version}")
    print()


def example_parse_sdf_string():
    """Parse a molecule from an SDF string."""
    print("=" * 60)
    print("Parsing SDF from String")
    print("=" * 60)

    # Simple methane molecule
    sdf_content = """\
methane


  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"""
    mol = sdfrust.parse_sdf_string(sdf_content)
    print(f"Parsed: {mol.name} with formula {mol.formula()}")
    print()


def example_create_molecule():
    """Create a molecule programmatically."""
    print("=" * 60)
    print("Creating Molecule Programmatically")
    print("=" * 60)

    # Create water molecule
    mol = sdfrust.Molecule("water")

    # Add atoms (index, element, x, y, z)
    mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
    mol.add_atom(sdfrust.Atom(1, "H", 0.9572, 0.0, 0.0))
    mol.add_atom(sdfrust.Atom(2, "H", -0.2400, 0.9266, 0.0))

    # Add bonds (atom1_index, atom2_index, bond_order)
    mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 2, sdfrust.BondOrder.single()))

    print(f"Created: {mol}")
    print(f"Formula: {mol.formula()}")
    print(f"Is empty: {mol.is_empty()}")
    print()


def example_atom_access():
    """Access and manipulate atoms."""
    print("=" * 60)
    print("Working with Atoms")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    # Access individual atom
    atom = mol.get_atom(0)
    print(f"First atom: {atom}")
    print(f"  Element: {atom.element}")
    print(f"  Coordinates: {atom.coords()}")
    print(f"  Formal charge: {atom.formal_charge}")
    print(f"  Is charged: {atom.is_charged()}")

    # Iterate over all atoms
    print(f"\nAll atoms ({mol.num_atoms} total):")
    for i, atom in enumerate(mol.atoms[:5]):  # First 5 atoms
        print(f"  {i}: {atom.element} at ({atom.x:.3f}, {atom.y:.3f}, {atom.z:.3f})")
    if mol.num_atoms > 5:
        print(f"  ... and {mol.num_atoms - 5} more")

    # Get atoms by element
    carbon_atoms = mol.atoms_by_element("C")
    oxygen_atoms = mol.atoms_by_element("O")
    print(f"\nCarbon atoms: {len(carbon_atoms)}")
    print(f"Oxygen atoms: {len(oxygen_atoms)}")
    print()


def example_bond_access():
    """Access and manipulate bonds."""
    print("=" * 60)
    print("Working with Bonds")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    # Access individual bond
    bond = mol.bonds[0]
    print(f"First bond: {bond}")
    print(f"  Connects atoms: {bond.atom1} - {bond.atom2}")
    print(f"  Order: {bond.order}")
    print(f"  Stereo: {bond.stereo}")

    # Bond order types
    print("\nBond order examples:")
    print(f"  Single: order = {sdfrust.BondOrder.single().order()}")
    print(f"  Double: order = {sdfrust.BondOrder.double().order()}")
    print(f"  Triple: order = {sdfrust.BondOrder.triple().order()}")
    print(f"  Aromatic: order = {sdfrust.BondOrder.aromatic().order()}")

    # Get bonds for a specific atom
    bonds_for_atom0 = mol.bonds_for_atom(0)
    print(f"\nBonds connected to atom 0: {len(bonds_for_atom0)}")

    # Get neighbors of an atom
    neighbors = mol.neighbors(0)
    print(f"Neighbors of atom 0: {neighbors}")

    # Filter bonds by type
    double_bonds = mol.bonds_by_order(sdfrust.BondOrder.double())
    print(f"\nDouble bonds in molecule: {len(double_bonds)}")

    # Bond type counts
    bond_counts = mol.bond_type_counts()
    print(f"Bond type distribution: {bond_counts}")
    print()


def example_molecular_descriptors():
    """Calculate molecular descriptors."""
    print("=" * 60)
    print("Molecular Descriptors")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    print(f"Molecule: {mol.name}")
    print(f"  Formula: {mol.formula()}")
    print(f"  Molecular weight: {mol.molecular_weight():.3f} g/mol")
    print(f"  Exact mass: {mol.exact_mass():.4f} g/mol")
    print(f"  Heavy atom count: {mol.heavy_atom_count()}")
    print(f"  Ring count: {mol.ring_count()}")
    print(f"  Rotatable bonds: {mol.rotatable_bond_count()}")
    print(f"  Total charge: {mol.total_charge()}")
    print(f"  Has charges: {mol.has_charges()}")
    print(f"  Has aromatic bonds: {mol.has_aromatic_bonds()}")
    print(f"  Total bond order: {mol.total_bond_order()}")

    # Element counts
    element_counts = mol.element_counts()
    print(f"  Element counts: {element_counts}")
    print()


def example_geometry():
    """Geometric operations on molecules."""
    print("=" * 60)
    print("Geometry Operations")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    # Calculate centroid
    centroid = mol.centroid()
    print(f"Centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")

    # Get first atom position
    atom0 = mol.get_atom(0)
    print(f"Atom 0 before translation: ({atom0.x:.3f}, {atom0.y:.3f}, {atom0.z:.3f})")

    # Translate molecule
    mol.translate(1.0, 2.0, 3.0)
    atom0 = mol.get_atom(0)
    print(f"Atom 0 after translate(1, 2, 3): ({atom0.x:.3f}, {atom0.y:.3f}, {atom0.z:.3f})")

    # Center molecule at origin
    mol.center()
    centroid = mol.centroid()
    print(f"Centroid after centering: ({centroid[0]:.6f}, {centroid[1]:.6f}, {centroid[2]:.6f})")

    # Distance between atoms
    atom1 = sdfrust.Atom(0, "C", 0.0, 0.0, 0.0)
    atom2 = sdfrust.Atom(1, "C", 3.0, 4.0, 0.0)
    print(f"\nDistance between (0,0,0) and (3,4,0): {atom1.distance_to(atom2):.3f}")
    print()


def example_properties():
    """Work with molecule properties from SDF data blocks."""
    print("=" * 60)
    print("Molecule Properties")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    # Get all properties
    props = mol.properties
    print(f"Properties in file ({len(props)} total):")
    for key, value in list(props.items())[:5]:  # First 5 properties
        # Truncate long values
        display_value = value[:50] + "..." if len(value) > 50 else value
        print(f"  {key}: {display_value}")

    # Get specific property
    if "PUBCHEM_COMPOUND_CID" in props:
        print(f"\nPubChem CID: {mol.get_property('PUBCHEM_COMPOUND_CID')}")

    # Set custom property
    mol.set_property("MY_PROPERTY", "custom value")
    print(f"Custom property: {mol.get_property('MY_PROPERTY')}")
    print()


def example_iterate_multi_molecule():
    """Iterate over multi-molecule SDF files."""
    print("=" * 60)
    print("Iterating Multi-Molecule Files")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")

    # Method 1: Load all molecules into a list
    molecules = sdfrust.parse_sdf_file_multi(aspirin_path)
    print(f"Loaded {len(molecules)} molecule(s) with parse_sdf_file_multi()")

    # Method 2: Memory-efficient streaming iterator
    print("\nUsing streaming iterator:")
    for i, mol in enumerate(sdfrust.iter_sdf_file(aspirin_path)):
        print(f"  Molecule {i}: {mol.name} ({mol.num_atoms} atoms)")

    # This is more memory-efficient for large files as it doesn't
    # load all molecules into memory at once
    print()


def example_write_sdf():
    """Write molecules to SDF files."""
    print("=" * 60)
    print("Writing SDF Files")
    print("=" * 60)

    # Create a simple molecule
    mol = sdfrust.Molecule("ethanol")
    # Ethanol: C2H5OH
    mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
    mol.add_atom(sdfrust.Atom(1, "C", 1.5, 0.0, 0.0))
    mol.add_atom(sdfrust.Atom(2, "O", 2.3, 1.0, 0.0))
    mol.add_atom(sdfrust.Atom(3, "H", -0.5, 0.9, 0.0))
    mol.add_atom(sdfrust.Atom(4, "H", -0.5, -0.9, 0.0))
    mol.add_atom(sdfrust.Atom(5, "H", 0.0, 0.0, 1.0))
    mol.add_atom(sdfrust.Atom(6, "H", 2.0, -0.9, 0.0))
    mol.add_atom(sdfrust.Atom(7, "H", 2.0, 0.9, 0.0))
    mol.add_atom(sdfrust.Atom(8, "H", 3.2, 0.8, 0.0))

    mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 2, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 3, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 4, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 5, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 6, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 7, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(2, 8, sdfrust.BondOrder.single()))

    mol.set_property("SMILES", "CCO")

    # Write to string (V2000 format)
    sdf_string = sdfrust.write_sdf_string(mol)
    print("SDF V2000 output (first 500 chars):")
    print(sdf_string[:500])
    print("...")

    # Write to string (V3000 format)
    sdf_v3000 = sdfrust.write_sdf_v3000_string(mol)
    print("\nSDF V3000 output (first 500 chars):")
    print(sdf_v3000[:500])
    print("...")

    # Write to file
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        temp_path = f.name
    sdfrust.write_sdf_file(mol, temp_path)
    print(f"\nWrote molecule to: {temp_path}")

    # Read back and verify
    mol_read = sdfrust.parse_sdf_file(temp_path)
    print(f"Read back: {mol_read.name} ({mol_read.formula()})")

    # Clean up
    os.unlink(temp_path)
    print()


def example_mol2_parsing():
    """Parse MOL2 (TRIPOS) format files."""
    print("=" * 60)
    print("Parsing MOL2 Files")
    print("=" * 60)

    benzene_path = os.path.join(TEST_DATA, "benzene.mol2")
    mol = sdfrust.parse_mol2_file(benzene_path)

    print(f"Molecule: {mol.name}")
    print(f"Atoms: {mol.num_atoms}")
    print(f"Bonds: {mol.num_bonds}")
    print(f"Formula: {mol.formula()}")

    # MOL2 can also be parsed from strings
    mol2_content = """\
@<TRIPOS>MOLECULE
water
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
"""
    mol = sdfrust.parse_mol2_string(mol2_content)
    print(f"\nParsed from string: {mol.name} ({mol.formula()})")
    print()


def example_v3000_format():
    """Work with SDF V3000 format."""
    print("=" * 60)
    print("SDF V3000 Format")
    print("=" * 60)

    v3000_path = os.path.join(TEST_DATA, "v3000_benzene.sdf")
    mol = sdfrust.parse_sdf_v3000_file(v3000_path)

    print(f"Molecule: {mol.name}")
    print(f"Format: {mol.format_version}")
    print(f"Atoms: {mol.num_atoms}")
    print(f"Formula: {mol.formula()}")

    # Check if molecule needs V3000 format
    print(f"Needs V3000: {mol.needs_v3000()}")

    # Auto-detect format when parsing
    mol_auto = sdfrust.parse_sdf_auto_file(v3000_path)
    print(f"\nAuto-detected format: {mol_auto.format_version}")
    print()


def example_numpy_integration():
    """NumPy integration for coordinate manipulation."""
    print("=" * 60)
    print("NumPy Integration")
    print("=" * 60)

    try:
        import numpy as np
    except ImportError:
        print("NumPy not available, skipping this example")
        print()
        return

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    mol = sdfrust.parse_sdf_file(aspirin_path)

    # Get coordinates as NumPy array
    coords = mol.get_coords_array()
    print(f"Coordinates shape: {coords.shape}")
    print(f"Coordinates dtype: {coords.dtype}")
    print(f"First 3 atoms:")
    print(coords[:3])

    # Get atomic numbers
    atomic_nums = mol.get_atomic_numbers()
    print(f"\nAtomic numbers: {atomic_nums[:10]}...")

    # Manipulate coordinates with NumPy
    # Example: Scale all coordinates by 2
    new_coords = coords * 2.0
    mol.set_coords_array(new_coords)

    # Verify the change
    atom0 = mol.get_atom(0)
    print(f"\nAfter scaling by 2:")
    print(f"  Original first coord: {coords[0]}")
    print(f"  New first coord: ({atom0.x:.3f}, {atom0.y:.3f}, {atom0.z:.3f})")
    print()


def example_ring_detection():
    """Detect rings and ring membership."""
    print("=" * 60)
    print("Ring Detection")
    print("=" * 60)

    # Load benzene (aromatic ring)
    benzene_path = os.path.join(TEST_DATA, "benzene.mol2")
    mol = sdfrust.parse_mol2_file(benzene_path)

    print(f"Molecule: {mol.name}")
    print(f"Ring count: {mol.ring_count()}")

    # Check which atoms are in rings
    print("\nAtoms in rings:")
    for i, atom in enumerate(mol.atoms):
        in_ring = mol.is_atom_in_ring(i)
        print(f"  Atom {i} ({atom.element}): {'in ring' if in_ring else 'not in ring'}")

    # Check which bonds are in rings
    print("\nBonds in rings:")
    for i, bond in enumerate(mol.bonds):
        in_ring = mol.is_bond_in_ring(i)
        print(f"  Bond {i} ({bond.atom1}-{bond.atom2}): {'in ring' if in_ring else 'not in ring'}")
    print()


def example_charged_atoms():
    """Work with charged atoms."""
    print("=" * 60)
    print("Charged Atoms")
    print("=" * 60)

    # Create a molecule with charged atoms
    mol = sdfrust.Molecule("ammonium")

    # NH4+ ion
    n_atom = sdfrust.Atom(0, "N", 0.0, 0.0, 0.0)
    n_atom.formal_charge = 1
    mol.add_atom(n_atom)

    mol.add_atom(sdfrust.Atom(1, "H", 1.0, 0.0, 0.0))
    mol.add_atom(sdfrust.Atom(2, "H", -0.33, 0.94, 0.0))
    mol.add_atom(sdfrust.Atom(3, "H", -0.33, -0.47, 0.82))
    mol.add_atom(sdfrust.Atom(4, "H", -0.33, -0.47, -0.82))

    mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 2, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 3, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 4, sdfrust.BondOrder.single()))

    print(f"Molecule: {mol.name}")
    print(f"Total charge: {mol.total_charge()}")
    print(f"Has charges: {mol.has_charges()}")

    # Check individual atoms
    for i, atom in enumerate(mol.atoms):
        print(f"  Atom {i} ({atom.element}): charge = {atom.formal_charge}")
    print()


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("  sdfrust Python Examples")
    print("=" * 60 + "\n")

    example_version()
    example_parse_sdf_file()
    example_parse_sdf_string()
    example_create_molecule()
    example_atom_access()
    example_bond_access()
    example_molecular_descriptors()
    example_geometry()
    example_properties()
    example_iterate_multi_molecule()
    example_write_sdf()
    example_mol2_parsing()
    example_v3000_format()
    example_numpy_integration()
    example_ring_detection()
    example_charged_atoms()

    print("=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
