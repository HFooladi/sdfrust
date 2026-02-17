#!/usr/bin/env python3
"""
Geometry Analysis Examples for sdfrust Python Bindings

This script demonstrates 3D geometry operations for structural analysis,
including distance matrices, RMSD comparison, rotation, transformation,
and NumPy-based coordinate analysis.

Requires the geometry feature:
    cd sdfrust-python
    maturin develop --features numpy,geometry
    python examples/geometry_analysis.py
"""

import math
import os

import sdfrust

# Path to test data shipped with the Rust crate
TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")

# Path to example-specific data (PubChem downloads)
EXAMPLE_DATA = os.path.join(os.path.dirname(__file__), "data")


def _has_geometry():
    """Check if the geometry feature is available."""
    return hasattr(sdfrust.Molecule, "distance_matrix")


def example_distance_matrix():
    """Compute pairwise distance matrix and find closest/farthest atom pairs."""
    print("=" * 60)
    print("Distance Matrix Analysis")
    print("=" * 60)

    # Use dopamine — small enough to display
    dopamine_path = os.path.join(EXAMPLE_DATA, "dopamine.sdf")
    mol = sdfrust.parse_sdf_file(dopamine_path)
    print(f"Molecule: {mol.name} ({mol.num_atoms} atoms)")

    matrix = mol.distance_matrix()
    n = mol.num_atoms

    # Print a compact distance matrix for heavy atoms only
    heavy_indices = [i for i, a in enumerate(mol.atoms) if a.element != "H"]
    print(f"\nDistance matrix (heavy atoms only, {len(heavy_indices)} atoms):")
    header = "      " + "".join(f"{mol.get_atom(j).element}{j:>3d}" for j in heavy_indices[:8])
    print(header)
    for i in heavy_indices[:8]:
        row = f"{mol.get_atom(i).element}{i:>3d}  "
        for j in heavy_indices[:8]:
            row += f"{matrix[i][j]:6.2f}"
        print(row)

    # Find closest non-bonded pair and farthest pair
    min_dist = float("inf")
    max_dist = 0.0
    min_pair = (0, 0)
    max_pair = (0, 0)
    for i in range(n):
        for j in range(i + 1, n):
            d = matrix[i][j]
            if d < min_dist:
                min_dist = d
                min_pair = (i, j)
            if d > max_dist:
                max_dist = d
                max_pair = (i, j)

    ai, aj = min_pair
    print(f"\nClosest pair:  atoms {ai} ({mol.get_atom(ai).element}) - "
          f"{aj} ({mol.get_atom(aj).element}) = {min_dist:.3f} A")
    ai, aj = max_pair
    print(f"Farthest pair: atoms {ai} ({mol.get_atom(ai).element}) - "
          f"{aj} ({mol.get_atom(aj).element}) = {max_dist:.3f} A")

    print()


def example_bond_length_analysis():
    """Extract actual bond lengths grouped by bond type."""
    print("=" * 60)
    print("Bond Length Analysis")
    print("=" * 60)

    dopamine_path = os.path.join(EXAMPLE_DATA, "dopamine.sdf")
    mol = sdfrust.parse_sdf_file(dopamine_path)
    matrix = mol.distance_matrix()

    print(f"Molecule: {mol.name}\n")

    # Group bonds by order and collect distances
    bond_lengths = {}
    for bond in mol.bonds:
        order_str = str(bond.order)
        d = matrix[bond.atom1][bond.atom2]
        if order_str not in bond_lengths:
            bond_lengths[order_str] = []
        bond_lengths[order_str].append((bond.atom1, bond.atom2, d))

    for order, bonds in sorted(bond_lengths.items()):
        distances = [d for _, _, d in bonds]
        avg = sum(distances) / len(distances)
        min_d = min(distances)
        max_d = max(distances)
        print(f"Bond order {order}: {len(bonds)} bonds")
        print(f"  Range: {min_d:.3f} - {max_d:.3f} A")
        print(f"  Mean:  {avg:.3f} A")
        # Show individual bonds
        for a1, a2, d in bonds[:5]:
            e1 = mol.get_atom(a1).element
            e2 = mol.get_atom(a2).element
            print(f"    {e1}{a1}-{e2}{a2}: {d:.3f} A")
        if len(bonds) > 5:
            print(f"    ... and {len(bonds) - 5} more")

    print()


def example_rmsd_comparison():
    """Compare RMSD between identical and translated molecules."""
    print("=" * 60)
    print("RMSD Comparison")
    print("=" * 60)

    ibu_path = os.path.join(EXAMPLE_DATA, "ibuprofen.sdf")
    mol1 = sdfrust.parse_sdf_file(ibu_path)
    mol2 = sdfrust.parse_sdf_file(ibu_path)

    # RMSD of identical molecules should be 0
    rmsd_same = mol1.rmsd_to(mol2)
    print(f"RMSD (identical copies): {rmsd_same:.6f} A")

    # Translate mol2 and check RMSD
    mol2.translate(1.0, 0.0, 0.0)
    rmsd_translated = mol1.rmsd_to(mol2)
    print(f"RMSD (after +1A x-shift): {rmsd_translated:.6f} A")

    # Translate in all directions
    mol3 = sdfrust.parse_sdf_file(ibu_path)
    mol3.translate(1.0, 2.0, 3.0)
    rmsd_3d = mol1.rmsd_to(mol3)
    expected = math.sqrt(1.0**2 + 2.0**2 + 3.0**2)
    print(f"RMSD (after +1,+2,+3 shift): {rmsd_3d:.6f} A (expected: {expected:.6f})")

    print()


def example_rotation():
    """Rotate a molecule and verify distances are preserved."""
    print("=" * 60)
    print("Rotation (Axis-Angle)")
    print("=" * 60)

    dopamine_path = os.path.join(EXAMPLE_DATA, "dopamine.sdf")
    mol = sdfrust.parse_sdf_file(dopamine_path)

    # Get distances before rotation
    dm_before = mol.distance_matrix()
    d01_before = dm_before[0][1]
    d02_before = dm_before[0][2]

    # Get centroid before
    centroid_before = mol.centroid()
    print(f"Centroid before: ({centroid_before[0]:.3f}, {centroid_before[1]:.3f}, {centroid_before[2]:.3f})")

    # Center molecule, then rotate 90 degrees around Z axis
    mol.center()
    mol.rotate([0.0, 0.0, 1.0], math.pi / 2)

    centroid_after = mol.centroid()
    print(f"Centroid after center+rotate: ({centroid_after[0]:.6f}, {centroid_after[1]:.6f}, {centroid_after[2]:.6f})")

    # Verify distances are preserved
    dm_after = mol.distance_matrix()
    d01_after = dm_after[0][1]
    d02_after = dm_after[0][2]

    print(f"\nDistance atom0-atom1: before={d01_before:.4f}, after={d01_after:.4f}")
    print(f"Distance atom0-atom2: before={d02_before:.4f}, after={d02_after:.4f}")

    diff01 = abs(d01_before - d01_after)
    diff02 = abs(d02_before - d02_after)
    print(f"Difference: {diff01:.2e}, {diff02:.2e}")
    print("Distances preserved: " + ("YES" if diff01 < 1e-6 and diff02 < 1e-6 else "NO"))

    print()


def example_apply_transform():
    """Apply combined rotation matrix + translation."""
    print("=" * 60)
    print("Combined Transform (Rotation Matrix + Translation)")
    print("=" * 60)

    dopamine_path = os.path.join(EXAMPLE_DATA, "dopamine.sdf")
    mol = sdfrust.parse_sdf_file(dopamine_path)

    atom0_before = mol.get_atom(0)
    print(f"Atom 0 before: ({atom0_before.x:.3f}, {atom0_before.y:.3f}, {atom0_before.z:.3f})")

    # Apply identity rotation with translation [10, 20, 30]
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    translation = [10.0, 20.0, 30.0]
    mol.apply_transform(identity, translation)

    atom0_after = mol.get_atom(0)
    print(f"Atom 0 after identity + [10,20,30]: ({atom0_after.x:.3f}, {atom0_after.y:.3f}, {atom0_after.z:.3f})")

    # Apply 180-degree rotation around Z axis (no translation)
    rot_180z = [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]
    mol2 = sdfrust.parse_sdf_file(dopamine_path)
    mol2.apply_rotation_matrix(rot_180z)

    atom0_rot = mol2.get_atom(0)
    print(f"Atom 0 after 180° Z rotation: ({atom0_rot.x:.3f}, {atom0_rot.y:.3f}, {atom0_rot.z:.3f})")

    # Verify: x and y should be negated, z unchanged
    orig = mol2_orig_atom0 = sdfrust.parse_sdf_file(dopamine_path).get_atom(0)
    print(f"  Expected: ({-orig.x:.3f}, {-orig.y:.3f}, {orig.z:.3f})")

    print()


def example_numpy_coordinate_analysis():
    """Use NumPy for bounding box, spatial extent, and radius of gyration."""
    print("=" * 60)
    print("NumPy Coordinate Analysis")
    print("=" * 60)

    try:
        import numpy as np
    except ImportError:
        print("NumPy not available, skipping this example")
        print()
        return

    chol_path = os.path.join(EXAMPLE_DATA, "cholesterol.sdf")
    mol = sdfrust.parse_sdf_file(chol_path)
    coords = mol.get_coords_array()
    print(f"Molecule: {mol.name} ({mol.num_atoms} atoms)")
    print(f"Coordinate array shape: {coords.shape}")

    # Bounding box
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    extent = max_coords - min_coords
    print(f"\nBounding box:")
    print(f"  Min: ({min_coords[0]:.3f}, {min_coords[1]:.3f}, {min_coords[2]:.3f})")
    print(f"  Max: ({max_coords[0]:.3f}, {max_coords[1]:.3f}, {max_coords[2]:.3f})")
    print(f"  Extent: ({extent[0]:.3f}, {extent[1]:.3f}, {extent[2]:.3f})")

    # Spatial extent (max distance from centroid)
    centroid = coords.mean(axis=0)
    distances_from_center = np.linalg.norm(coords - centroid, axis=1)
    max_extent = distances_from_center.max()
    print(f"\nSpatial extent from centroid: {max_extent:.3f} A")
    print(f"Mean distance from centroid:  {distances_from_center.mean():.3f} A")

    # Radius of gyration
    rg = np.sqrt(np.mean(distances_from_center**2))
    print(f"Radius of gyration: {rg:.3f} A")

    # Per-element spatial analysis
    print("\nPer-element centroid offsets:")
    atomic_nums = mol.get_atomic_numbers()
    for elem_str, _ in sorted(mol.element_counts().items()):
        elem_atoms = mol.atoms_by_element(elem_str)
        if len(elem_atoms) > 0:
            elem_coords = np.array([[a.x, a.y, a.z] for a in elem_atoms])
            elem_centroid = elem_coords.mean(axis=0)
            offset = np.linalg.norm(elem_centroid - centroid)
            print(f"  {elem_str:>2s} ({len(elem_atoms):3d} atoms): offset from centroid = {offset:.3f} A")

    print()


def example_compare_conformers():
    """Simulate conformer comparison: perturb, center, compare RMSD."""
    print("=" * 60)
    print("Conformer Comparison (Perturb + Center + RMSD)")
    print("=" * 60)

    try:
        import numpy as np
    except ImportError:
        print("NumPy not available, skipping this example")
        print()
        return

    ibu_path = os.path.join(EXAMPLE_DATA, "ibuprofen.sdf")
    mol1 = sdfrust.parse_sdf_file(ibu_path)
    mol2 = sdfrust.parse_sdf_file(ibu_path)

    # Add small random perturbation to mol2 coordinates
    np.random.seed(42)
    coords = mol2.get_coords_array()
    perturbation = np.random.normal(0, 0.1, coords.shape)  # 0.1 A std dev
    mol2.set_coords_array(coords + perturbation)

    # RMSD before centering
    rmsd_before = mol1.rmsd_to(mol2)
    print(f"RMSD (before centering): {rmsd_before:.4f} A")

    # Center both molecules
    mol1.center()
    mol2.center()

    # RMSD after centering
    rmsd_after = mol1.rmsd_to(mol2)
    print(f"RMSD (after centering):  {rmsd_after:.4f} A")

    # The RMSD should be similar (centering removes translational component
    # but the perturbation is the dominant effect here)
    print(f"Perturbation std dev:    0.100 A")
    print(f"Expected RMSD ~ 0.1 * sqrt(3) = {0.1 * math.sqrt(3):.4f} A")

    print()


def main():
    """Run all geometry analysis examples."""
    print("\n" + "=" * 60)
    print("  sdfrust Geometry Analysis Examples")
    print("=" * 60 + "\n")

    if not _has_geometry():
        print("NOTICE: Geometry feature not available.")
        print("Rebuild with geometry feature enabled:")
        print('  maturin develop --features "numpy,geometry"')
        print()
        print("Skipping geometry-specific examples.")
        print("Running NumPy-only examples...\n")
        example_numpy_coordinate_analysis()
        example_compare_conformers()
        print("=" * 60)
        print("Geometry examples completed (partial — no geometry feature)!")
        print("=" * 60)
        return

    example_distance_matrix()
    example_bond_length_analysis()
    example_rmsd_comparison()
    example_rotation()
    example_apply_transform()
    example_numpy_coordinate_analysis()
    example_compare_conformers()

    print("=" * 60)
    print("All geometry analysis examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
