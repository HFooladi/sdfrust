#!/usr/bin/env python3
"""
Format Conversion Examples for sdfrust Python Bindings

This script demonstrates multi-format detection, parsing, and conversion
capabilities, including XYZ parsing, MOL2 writing, batch conversions,
and round-trip verification.

Run this script from the sdfrust-python directory after building:
    cd sdfrust-python
    maturin develop --features numpy
    python examples/format_conversion.py
"""

import os
import tempfile

import sdfrust

# Path to test data shipped with the Rust crate
TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")

# Path to example-specific data (PubChem downloads)
EXAMPLE_DATA = os.path.join(os.path.dirname(__file__), "data")


def example_detect_format():
    """Detect file format from content."""
    print("=" * 60)
    print("Detecting File Formats")
    print("=" * 60)

    # Read SDF file content and detect format
    sdf_path = os.path.join(EXAMPLE_DATA, "ibuprofen.sdf")
    with open(sdf_path) as f:
        sdf_content = f.read()
    fmt = sdfrust.detect_format(sdf_content)
    print(f"Ibuprofen SDF detected as: {fmt}")

    # Read MOL2 file content and detect format
    mol2_path = os.path.join(TEST_DATA, "benzene.mol2")
    with open(mol2_path) as f:
        mol2_content = f.read()
    fmt = sdfrust.detect_format(mol2_content)
    print(f"Benzene MOL2 detected as:  {fmt}")

    # Read V3000 file content and detect format
    v3000_path = os.path.join(TEST_DATA, "v3000_benzene.sdf")
    with open(v3000_path) as f:
        v3000_content = f.read()
    fmt = sdfrust.detect_format(v3000_content)
    print(f"V3000 benzene detected as: {fmt}")

    print()


def example_parse_auto_file():
    """Parse any supported format with a single function."""
    print("=" * 60)
    print("Auto-Detection Parsing (parse_auto_file)")
    print("=" * 60)

    # Parse SDF file — auto-detected
    sdf_path = os.path.join(EXAMPLE_DATA, "ibuprofen.sdf")
    mol = sdfrust.parse_auto_file(sdf_path)
    print(f"SDF:  {mol.name:>12s} | {mol.num_atoms:3d} atoms | {mol.formula()}")

    # Parse MOL2 file — auto-detected
    mol2_path = os.path.join(TEST_DATA, "benzene.mol2")
    mol = sdfrust.parse_auto_file(mol2_path)
    print(f"MOL2: {mol.name:>12s} | {mol.num_atoms:3d} atoms | {mol.formula()}")

    # parse_auto_file handles V3000 as well
    v3000_path = os.path.join(TEST_DATA, "v3000_benzene.sdf")
    mol = sdfrust.parse_auto_file(v3000_path)
    print(f"V3K:  {mol.name:>12s} | {mol.num_atoms:3d} atoms | {mol.formula()}")

    print()


def example_parse_xyz():
    """Parse XYZ format files — single and multi-molecule."""
    print("=" * 60)
    print("XYZ Format Parsing")
    print("=" * 60)

    # Single molecule XYZ
    water_path = os.path.join(TEST_DATA, "water.xyz")
    mol = sdfrust.parse_xyz_file(water_path)
    print(f"Single XYZ: {mol.name}")
    print(f"  Atoms: {mol.num_atoms}")
    print(f"  Bonds: {mol.num_bonds}  (XYZ has no bond info)")
    print(f"  Formula: {mol.formula()}")

    # Multi-molecule XYZ (bulk load)
    multi_path = os.path.join(TEST_DATA, "multi.xyz")
    molecules = sdfrust.parse_xyz_file_multi(multi_path)
    print(f"\nMulti XYZ: loaded {len(molecules)} molecules")
    for i, mol in enumerate(molecules):
        print(f"  [{i}] {mol.name}: {mol.num_atoms} atoms, {mol.formula()}")

    # Streaming iterator for XYZ
    print(f"\nStreaming XYZ with iter_xyz_file():")
    for i, mol in enumerate(sdfrust.iter_xyz_file(multi_path)):
        centroid = mol.centroid()
        cx, cy, cz = centroid if centroid else (0, 0, 0)
        print(f"  [{i}] {mol.name}: centroid = ({cx:.3f}, {cy:.3f}, {cz:.3f})")

    print()


def example_sdf_to_mol2():
    """Convert an SDF molecule to MOL2 format with round-trip verification."""
    print("=" * 60)
    print("SDF → MOL2 Conversion (Round-Trip)")
    print("=" * 60)

    # Parse ibuprofen from SDF
    sdf_path = os.path.join(EXAMPLE_DATA, "ibuprofen.sdf")
    mol_original = sdfrust.parse_sdf_file(sdf_path)
    print(f"Original (SDF):  {mol_original.name} | "
          f"{mol_original.num_atoms} atoms | {mol_original.num_bonds} bonds")

    # Write to MOL2 string
    mol2_string = sdfrust.write_mol2_string(mol_original)
    print(f"MOL2 output: {len(mol2_string)} characters")

    # Write to MOL2 file and read back
    with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as f:
        temp_mol2 = f.name
    try:
        sdfrust.write_mol2_file(mol_original, temp_mol2)
        mol_roundtrip = sdfrust.parse_mol2_file(temp_mol2)
        print(f"Round-trip (MOL2): {mol_roundtrip.name} | "
              f"{mol_roundtrip.num_atoms} atoms | {mol_roundtrip.num_bonds} bonds")

        # Verify atom count preserved
        assert mol_original.num_atoms == mol_roundtrip.num_atoms, "Atom count mismatch!"
        assert mol_original.num_bonds == mol_roundtrip.num_bonds, "Bond count mismatch!"
        print("Round-trip verification: PASSED (atom/bond counts match)")
    finally:
        os.unlink(temp_mol2)

    print()


def example_mol2_to_sdf():
    """Convert a MOL2 molecule to SDF format with round-trip verification."""
    print("=" * 60)
    print("MOL2 → SDF Conversion (Round-Trip)")
    print("=" * 60)

    # Parse benzene from MOL2
    mol2_path = os.path.join(TEST_DATA, "benzene.mol2")
    mol_original = sdfrust.parse_mol2_file(mol2_path)
    print(f"Original (MOL2): {mol_original.name} | "
          f"{mol_original.num_atoms} atoms | {mol_original.num_bonds} bonds")

    # Write to SDF string
    sdf_string = sdfrust.write_sdf_string(mol_original)
    print(f"SDF output: {len(sdf_string)} characters")

    # Write to SDF file and read back
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        temp_sdf = f.name
    try:
        sdfrust.write_sdf_file(mol_original, temp_sdf)
        mol_roundtrip = sdfrust.parse_sdf_file(temp_sdf)
        print(f"Round-trip (SDF):  {mol_roundtrip.name} | "
              f"{mol_roundtrip.num_atoms} atoms | {mol_roundtrip.num_bonds} bonds")

        assert mol_original.num_atoms == mol_roundtrip.num_atoms, "Atom count mismatch!"
        assert mol_original.num_bonds == mol_roundtrip.num_bonds, "Bond count mismatch!"
        print("Round-trip verification: PASSED (atom/bond counts match)")
    finally:
        os.unlink(temp_sdf)

    print()


def example_batch_format_conversion():
    """Convert multiple molecules between formats in batch."""
    print("=" * 60)
    print("Batch Format Conversion")
    print("=" * 60)

    # Load multiple molecules from different SDF files
    molecules = []
    for name in ["ibuprofen.sdf", "dopamine.sdf", "cholesterol.sdf"]:
        path = os.path.join(EXAMPLE_DATA, name)
        mol = sdfrust.parse_sdf_file(path)
        molecules.append(mol)
        print(f"Loaded: {mol.name:>12s} | {mol.formula()}")

    # Write all molecules to a single multi-molecule SDF file
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        temp_sdf = f.name
    try:
        sdfrust.write_sdf_file_multi(molecules, temp_sdf)
        print(f"\nWrote {len(molecules)} molecules to multi-SDF: {os.path.basename(temp_sdf)}")

        # Read back and verify
        mols_back = sdfrust.parse_sdf_file_multi(temp_sdf)
        print(f"Read back {len(mols_back)} molecules from multi-SDF")
        for mol in mols_back:
            print(f"  {mol.name}: {mol.num_atoms} atoms")
    finally:
        os.unlink(temp_sdf)

    # Write all molecules to a single multi-molecule MOL2 file
    with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as f:
        temp_mol2 = f.name
    try:
        sdfrust.write_mol2_file_multi(molecules, temp_mol2)
        print(f"\nWrote {len(molecules)} molecules to multi-MOL2: {os.path.basename(temp_mol2)}")

        # Read back and verify
        mols_back = sdfrust.parse_mol2_file_multi(temp_mol2)
        print(f"Read back {len(mols_back)} molecules from multi-MOL2")
        for mol in mols_back:
            print(f"  {mol.name}: {mol.num_atoms} atoms")
    finally:
        os.unlink(temp_mol2)

    print()


def example_gzip_support():
    """Check gzip runtime feature availability."""
    print("=" * 60)
    print("Gzip Support")
    print("=" * 60)

    if sdfrust.gzip_enabled():
        print("Gzip support: ENABLED")
        print("  You can parse .sdf.gz, .mol2.gz, and .xyz.gz files directly:")
        print('  mol = sdfrust.parse_sdf_file("molecule.sdf.gz")')
    else:
        print("Gzip support: DISABLED")
        print("  Rebuild with gzip feature to enable:")
        print('  maturin develop --features "numpy,gzip"')

    print()


def main():
    """Run all format conversion examples."""
    print("\n" + "=" * 60)
    print("  sdfrust Format Conversion Examples")
    print("=" * 60 + "\n")

    example_detect_format()
    example_parse_auto_file()
    example_parse_xyz()
    example_sdf_to_mol2()
    example_mol2_to_sdf()
    example_batch_format_conversion()
    example_gzip_support()

    print("=" * 60)
    print("All format conversion examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
