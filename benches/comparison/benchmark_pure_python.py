#!/usr/bin/env python3
"""
Pure Python SDF parser benchmark (no external dependencies except psutil).

This provides a baseline comparison for what pure Python parsing looks like
without optimized C/C++ libraries.

Usage:
    python benchmark_pure_python.py <sdf_file> [--output results.json]
"""

import argparse
import json
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

try:
    import psutil
except ImportError:
    print("Error: psutil not installed. Run: pip install psutil")
    sys.exit(1)


@dataclass
class Atom:
    """Represents an atom with coordinates and element."""
    index: int
    element: str
    x: float
    y: float
    z: float
    charge: int = 0


@dataclass
class Bond:
    """Represents a bond between two atoms."""
    atom1: int
    atom2: int
    order: int
    stereo: int = 0


@dataclass
class Molecule:
    """Container for molecule data."""
    name: str = ""
    atoms: list = field(default_factory=list)
    bonds: list = field(default_factory=list)
    properties: dict = field(default_factory=dict)


def parse_sdf_molecule(lines: list[str]) -> Optional[Molecule]:
    """Parse a single molecule from SDF lines."""
    if not lines:
        return None

    mol = Molecule()

    # Line 1: Molecule name
    mol.name = lines[0].strip() if lines else ""

    # Lines 2-3: Comment/program info (skip)
    # Line 4: Counts line
    if len(lines) < 4:
        return None

    counts_line = lines[3]
    try:
        atom_count = int(counts_line[0:3].strip())
        bond_count = int(counts_line[3:6].strip())
    except (ValueError, IndexError):
        return None

    # Parse atom block
    line_idx = 4
    for i in range(atom_count):
        if line_idx >= len(lines):
            break
        line = lines[line_idx]
        try:
            x = float(line[0:10].strip())
            y = float(line[10:20].strip())
            z = float(line[20:30].strip())
            element = line[31:34].strip()

            # Parse charge (positions 36-38)
            charge_code = 0
            if len(line) > 38:
                try:
                    charge_code = int(line[36:39].strip())
                except ValueError:
                    pass

            # Convert charge code to actual charge
            charge_map = {0: 0, 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3}
            charge = charge_map.get(charge_code, 0)

            mol.atoms.append(Atom(i, element, x, y, z, charge))
        except (ValueError, IndexError):
            pass
        line_idx += 1

    # Parse bond block
    for i in range(bond_count):
        if line_idx >= len(lines):
            break
        line = lines[line_idx]
        try:
            atom1 = int(line[0:3].strip()) - 1  # Convert to 0-based
            atom2 = int(line[3:6].strip()) - 1
            order = int(line[6:9].strip())
            stereo = 0
            if len(line) > 11:
                try:
                    stereo = int(line[9:12].strip())
                except ValueError:
                    pass
            mol.bonds.append(Bond(atom1, atom2, order, stereo))
        except (ValueError, IndexError):
            pass
        line_idx += 1

    # Parse M lines and properties
    while line_idx < len(lines):
        line = lines[line_idx]
        if line.startswith("M  END"):
            line_idx += 1
            break
        line_idx += 1

    # Parse data block (properties)
    while line_idx < len(lines):
        line = lines[line_idx]
        if line.startswith("> <"):
            # Extract property name
            match = re.match(r"> <([^>]+)>", line)
            if match:
                prop_name = match.group(1)
                line_idx += 1
                values = []
                while line_idx < len(lines) and lines[line_idx].strip():
                    values.append(lines[line_idx].strip())
                    line_idx += 1
                mol.properties[prop_name] = "\n".join(values)
        line_idx += 1

    return mol


def iter_sdf_file(filepath: str) -> Iterator[Molecule]:
    """Iterate over molecules in an SDF file."""
    current_lines = []

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("$$$$"):
                mol = parse_sdf_molecule(current_lines)
                if mol:
                    yield mol
                current_lines = []
            else:
                current_lines.append(line.rstrip("\n\r"))

    # Handle last molecule if no trailing $$$$
    if current_lines:
        mol = parse_sdf_molecule(current_lines)
        if mol:
            yield mol


def get_memory_mb():
    """Get current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def benchmark_sdf_parse(sdf_file: str, iterations: int = 5):
    """Benchmark SDF parsing with pure Python parser."""
    results = {
        "tool": "pure_python",
        "file": str(sdf_file),
        "iterations": iterations,
        "runs": [],
    }

    for i in range(iterations):
        # Force garbage collection before measurement
        import gc
        gc.collect()

        memory_before = get_memory_mb()
        peak_memory = memory_before

        start_time = time.perf_counter()
        first_mol_time = None
        mol_count = 0

        for mol in iter_sdf_file(sdf_file):
            if first_mol_time is None:
                first_mol_time = time.perf_counter() - start_time
            mol_count += 1

            # Sample memory periodically
            if mol_count % 1000 == 0:
                current_memory = get_memory_mb()
                peak_memory = max(peak_memory, current_memory)

        end_time = time.perf_counter()
        elapsed = end_time - start_time

        # Final memory measurement
        memory_after = get_memory_mb()
        peak_memory = max(peak_memory, memory_after)

        run_result = {
            "run": i + 1,
            "molecules": mol_count,
            "elapsed_seconds": elapsed,
            "molecules_per_second": mol_count / elapsed if elapsed > 0 else 0,
            "first_molecule_latency_ms": (first_mol_time * 1000) if first_mol_time else None,
            "memory_before_mb": memory_before,
            "memory_after_mb": memory_after,
            "peak_memory_mb": peak_memory,
        }
        results["runs"].append(run_result)

        print(f"Run {i + 1}/{iterations}: {mol_count} molecules in {elapsed:.3f}s "
              f"({run_result['molecules_per_second']:.0f} mol/s)")

    # Calculate averages
    if results["runs"]:
        results["average"] = {
            "molecules_per_second": sum(r["molecules_per_second"] for r in results["runs"]) / len(results["runs"]),
            "first_molecule_latency_ms": sum(r["first_molecule_latency_ms"] or 0 for r in results["runs"]) / len(results["runs"]),
            "peak_memory_mb": max(r["peak_memory_mb"] for r in results["runs"]),
        }

    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark pure Python SDF parsing")
    parser.add_argument("sdf_file", help="Path to SDF file to benchmark")
    parser.add_argument("--iterations", "-n", type=int, default=5,
                        help="Number of benchmark iterations (default: 5)")
    parser.add_argument("--output", "-o", help="Output JSON file for results")
    args = parser.parse_args()

    sdf_path = Path(args.sdf_file)
    if not sdf_path.exists():
        print(f"Error: File not found: {sdf_path}")
        sys.exit(1)

    print(f"Benchmarking Pure Python parser on: {sdf_path}")
    print(f"File size: {sdf_path.stat().st_size / (1024 * 1024):.2f} MB")
    print("-" * 60)

    results = benchmark_sdf_parse(str(sdf_path), args.iterations)

    print("-" * 60)
    print(f"Average: {results['average']['molecules_per_second']:.0f} molecules/second")
    print(f"First molecule latency: {results['average']['first_molecule_latency_ms']:.2f} ms")
    print(f"Peak memory: {results['average']['peak_memory_mb']:.1f} MB")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to: {args.output}")

    return results


if __name__ == "__main__":
    main()
