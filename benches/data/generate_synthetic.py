#!/usr/bin/env python3
"""
Generate synthetic SDF files for benchmarking.

This script generates valid SDF files with configurable numbers of molecules.
Each molecule is a simple linear carbon chain to ensure consistent structure
across benchmarks.

Usage:
    python generate_synthetic.py -n 10000 -o synthetic_10000.sdf
"""

import argparse
import random
from pathlib import Path


def generate_chain_molecule(name: str, num_atoms: int = 5, include_properties: bool = True) -> str:
    """Generate a linear carbon chain molecule in SDF format."""
    num_bonds = num_atoms - 1 if num_atoms > 1 else 0

    lines = []

    # Header block (3 lines)
    lines.append(name)
    lines.append("  sdfrust_benchmark")
    lines.append("")

    # Counts line
    lines.append(f"{num_atoms:3d}{num_bonds:3d}  0  0  0  0  0  0  0  0999 V2000")

    # Atom block - linear chain
    for i in range(num_atoms):
        x = i * 1.54  # C-C bond length approximation
        y = 0.0
        z = 0.0
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0")

    # Bond block - single bonds connecting consecutive atoms
    for i in range(num_bonds):
        atom1 = i + 1  # 1-based indexing
        atom2 = i + 2
        lines.append(f"{atom1:3d}{atom2:3d}  1  0  0  0  0")

    lines.append("M  END")

    # Properties (optional)
    if include_properties:
        lines.append("> <BENCHMARK_ID>")
        lines.append(name)
        lines.append("")

        lines.append("> <ATOM_COUNT>")
        lines.append(str(num_atoms))
        lines.append("")

    # Record terminator
    lines.append("$$$$")

    return "\n".join(lines) + "\n"


def generate_varied_molecule(mol_id: int) -> str:
    """Generate a molecule with varied size (5-30 atoms)."""
    num_atoms = random.randint(5, 30)
    name = f"mol_{mol_id:08d}"
    return generate_chain_molecule(name, num_atoms, include_properties=True)


def generate_sdf_file(output_path: str, num_molecules: int, varied: bool = False):
    """Generate an SDF file with the specified number of molecules."""
    print(f"Generating {num_molecules} molecules...")

    with open(output_path, 'w') as f:
        for i in range(num_molecules):
            if varied:
                mol = generate_varied_molecule(i)
            else:
                # Default: uniform 10-atom molecules
                name = f"mol_{i:08d}"
                mol = generate_chain_molecule(name, num_atoms=10)

            f.write(mol)

            if (i + 1) % 10000 == 0:
                print(f"  Generated {i + 1} molecules...")

    # Report file stats
    path = Path(output_path)
    size_mb = path.stat().st_size / (1024 * 1024)
    print(f"Done! Generated {output_path}")
    print(f"  - Molecules: {num_molecules}")
    print(f"  - File size: {size_mb:.2f} MB")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic SDF files for benchmarking")
    parser.add_argument("-n", "--num-molecules", type=int, default=10000,
                        help="Number of molecules to generate (default: 10000)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output SDF file path")
    parser.add_argument("--varied", action="store_true",
                        help="Generate molecules with varied sizes (5-30 atoms)")
    parser.add_argument("--atoms", type=int, default=10,
                        help="Number of atoms per molecule (default: 10, ignored if --varied)")
    args = parser.parse_args()

    generate_sdf_file(args.output, args.num_molecules, varied=args.varied)


if __name__ == "__main__":
    main()
