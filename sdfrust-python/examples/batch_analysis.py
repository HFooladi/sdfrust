#!/usr/bin/env python3
"""
Batch Analysis Examples for sdfrust Python Bindings

This script demonstrates realistic compound library processing — filtering,
sorting, and comparing molecules by descriptors using a multi-molecule SDF
file containing real drug molecules from PubChem.

Run this script from the sdfrust-python directory after building:
    cd sdfrust-python
    maturin develop --features numpy
    python examples/batch_analysis.py
"""

import os
import tempfile

import sdfrust

# Path to test data shipped with the Rust crate
TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")

# Path to example-specific data (PubChem downloads)
EXAMPLE_DATA = os.path.join(os.path.dirname(__file__), "data")


def example_load_drug_library():
    """Load a multi-molecule SDF file and print a summary table."""
    print("=" * 60)
    print("Loading Drug Library")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)
    print(f"Loaded {len(molecules)} molecules from drug_library.sdf\n")

    # Print summary table
    header = f"{'Name':>15s} {'Formula':>12s} {'MW':>8s} {'Atoms':>6s} {'Bonds':>6s} {'Rings':>6s} {'RotBonds':>9s}"
    print(header)
    print("-" * len(header))
    for mol in molecules:
        mw = mol.molecular_weight()
        mw_str = f"{mw:.1f}" if mw is not None else "N/A"
        print(f"{mol.name:>15s} {mol.formula():>12s} {mw_str:>8s} "
              f"{mol.num_atoms:>6d} {mol.num_bonds:>6d} "
              f"{mol.ring_count():>6d} {mol.rotatable_bond_count():>9d}")

    print()
    return molecules


def example_streaming_iteration():
    """Process molecules one at a time with memory-efficient iteration."""
    print("=" * 60)
    print("Streaming Iteration (iter_sdf_file)")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")

    total_atoms = 0
    total_bonds = 0
    count = 0

    print("Processing molecules one at a time:")
    for mol in sdfrust.iter_sdf_file(library_path):
        total_atoms += mol.num_atoms
        total_bonds += mol.num_bonds
        count += 1
        mw = mol.molecular_weight()
        mw_str = f"{mw:.1f}" if mw is not None else "N/A"
        print(f"  [{count}] {mol.name}: {mol.num_atoms} atoms, MW = {mw_str}")

    print(f"\nSummary: {count} molecules, {total_atoms} total atoms, {total_bonds} total bonds")
    print(f"Average: {total_atoms / count:.1f} atoms/molecule, "
          f"{total_bonds / count:.1f} bonds/molecule")
    print()


def example_filter_by_molecular_weight():
    """Filter molecules by molecular weight range and write subset to file."""
    print("=" * 60)
    print("Filtering by Molecular Weight (150-300 Da)")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    # Filter to MW range 150-300
    mw_min, mw_max = 150.0, 300.0
    filtered = []
    for mol in molecules:
        mw = mol.molecular_weight()
        if mw is not None and mw_min <= mw <= mw_max:
            filtered.append(mol)

    print(f"Total molecules: {len(molecules)}")
    print(f"In MW range [{mw_min:.0f}, {mw_max:.0f}]: {len(filtered)}")
    for mol in filtered:
        mw = mol.molecular_weight()
        print(f"  {mol.name}: MW = {mw:.2f}")

    # Write filtered subset to a new file
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        temp_path = f.name
    try:
        sdfrust.write_sdf_file_multi(filtered, temp_path)
        print(f"\nWrote {len(filtered)} molecules to {os.path.basename(temp_path)}")

        # Verify by reading back
        mols_back = sdfrust.parse_sdf_file_multi(temp_path)
        print(f"Verified: read back {len(mols_back)} molecules")
    finally:
        os.unlink(temp_path)

    print()


def example_sort_by_descriptor():
    """Sort molecules by molecular weight and heavy atom count."""
    print("=" * 60)
    print("Sorting by Descriptors")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    # Sort by molecular weight (ascending)
    sorted_by_mw = sorted(molecules, key=lambda m: m.molecular_weight() or 0)
    print("Ranked by Molecular Weight (ascending):")
    for rank, mol in enumerate(sorted_by_mw, 1):
        mw = mol.molecular_weight()
        mw_str = f"{mw:.2f}" if mw is not None else "N/A"
        print(f"  {rank}. {mol.name}: MW = {mw_str}")

    # Sort by heavy atom count (descending)
    sorted_by_heavy = sorted(molecules, key=lambda m: m.heavy_atom_count(), reverse=True)
    print("\nRanked by Heavy Atom Count (descending):")
    for rank, mol in enumerate(sorted_by_heavy, 1):
        print(f"  {rank}. {mol.name}: {mol.heavy_atom_count()} heavy atoms")

    print()


def example_compare_element_compositions():
    """Compare elemental compositions across molecules."""
    print("=" * 60)
    print("Element Composition Comparison")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    # Collect all elements present in the library
    all_elements = set()
    for mol in molecules:
        all_elements.update(mol.element_counts().keys())

    # Sort elements: C, H first, then alphabetically
    priority = {"C": 0, "H": 1}
    sorted_elements = sorted(all_elements, key=lambda e: (priority.get(e, 2), e))

    # Print composition table
    header = f"{'Molecule':>15s}" + "".join(f"{e:>5s}" for e in sorted_elements)
    print(header)
    print("-" * len(header))
    for mol in molecules:
        counts = mol.element_counts()
        row = f"{mol.name:>15s}"
        for elem in sorted_elements:
            count = counts.get(elem, 0)
            row += f"{count:>5d}"
        print(row)

    print()


def example_compare_bond_profiles():
    """Compare bond type distributions across molecules."""
    print("=" * 60)
    print("Bond Profile Comparison")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    # Collect all bond types
    all_types = set()
    for mol in molecules:
        all_types.update(mol.bond_type_counts().keys())
    sorted_types = sorted(all_types)

    # Print bond profile table
    header = f"{'Molecule':>15s}" + "".join(f"{bt:>10s}" for bt in sorted_types) + f"{'Aromatic?':>10s}"
    print(header)
    print("-" * len(header))
    for mol in molecules:
        counts = mol.bond_type_counts()
        row = f"{mol.name:>15s}"
        for bt in sorted_types:
            row += f"{counts.get(bt, 0):>10d}"
        row += f"{'Yes':>10s}" if mol.has_aromatic_bonds() else f"{'No':>10s}"
        print(row)

    print()


def example_property_analysis():
    """Access PubChem SDF properties and apply Lipinski's Rule of Five."""
    print("=" * 60)
    print("Property Analysis (Lipinski's Rule of Five)")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    print("Lipinski criteria: MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5\n")

    for mol in molecules:
        mw = mol.molecular_weight() or 0.0

        # Use element counts as proxy for H-bond donors/acceptors
        # HBA ~ count of N + O atoms
        # HBD ~ approximated from NH + OH groups (simplified: N + O count)
        elem_counts = mol.element_counts()
        n_count = elem_counts.get("N", 0)
        o_count = elem_counts.get("O", 0)
        hba = n_count + o_count  # Simplified HBA estimate

        # Check Lipinski criteria (simplified — no LogP without force field)
        mw_ok = mw <= 500
        hba_ok = hba <= 10

        violations = 0
        if not mw_ok:
            violations += 1
        if not hba_ok:
            violations += 1

        status = "PASS" if violations == 0 else f"WARN ({violations} violation{'s' if violations > 1 else ''})"
        print(f"  {mol.name:>15s}: MW={mw:7.1f}  HBA(N+O)={hba:2d}  → {status}")

    print()


def example_find_extremes():
    """Find molecules with extreme descriptor values."""
    print("=" * 60)
    print("Finding Extremes in the Library")
    print("=" * 60)

    library_path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    molecules = sdfrust.parse_sdf_file_multi(library_path)

    # Most atoms
    most_atoms = max(molecules, key=lambda m: m.num_atoms)
    print(f"Most atoms:      {most_atoms.name} ({most_atoms.num_atoms} atoms)")

    # Fewest atoms
    fewest_atoms = min(molecules, key=lambda m: m.num_atoms)
    print(f"Fewest atoms:    {fewest_atoms.name} ({fewest_atoms.num_atoms} atoms)")

    # Highest MW
    highest_mw = max(molecules, key=lambda m: m.molecular_weight() or 0)
    print(f"Highest MW:      {highest_mw.name} ({highest_mw.molecular_weight():.1f} Da)")

    # Lowest MW
    lowest_mw = min(molecules, key=lambda m: m.molecular_weight() or float('inf'))
    print(f"Lowest MW:       {lowest_mw.name} ({lowest_mw.molecular_weight():.1f} Da)")

    # Most rings
    most_rings = max(molecules, key=lambda m: m.ring_count())
    print(f"Most rings:      {most_rings.name} ({most_rings.ring_count()} rings)")

    # Most rotatable bonds
    most_rotatable = max(molecules, key=lambda m: m.rotatable_bond_count())
    print(f"Most rot. bonds: {most_rotatable.name} ({most_rotatable.rotatable_bond_count()} rotatable bonds)")

    # Most bonds
    most_bonds = max(molecules, key=lambda m: m.num_bonds)
    print(f"Most bonds:      {most_bonds.name} ({most_bonds.num_bonds} bonds)")

    print()


def main():
    """Run all batch analysis examples."""
    print("\n" + "=" * 60)
    print("  sdfrust Batch Analysis Examples")
    print("=" * 60 + "\n")

    example_load_drug_library()
    example_streaming_iteration()
    example_filter_by_molecular_weight()
    example_sort_by_descriptor()
    example_compare_element_compositions()
    example_compare_bond_profiles()
    example_property_analysis()
    example_find_extremes()

    print("=" * 60)
    print("All batch analysis examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
