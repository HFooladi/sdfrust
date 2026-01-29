#!/usr/bin/env python3
"""
Benchmark RDKit SDF parsing performance.

This script measures:
- Parse throughput (molecules/second)
- Peak memory usage
- First molecule latency

Usage:
    python benchmark_rdkit.py <sdf_file> [--output results.json]
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

try:
    import psutil
except ImportError:
    print("Error: psutil not installed. Run: pip install psutil")
    sys.exit(1)

try:
    from rdkit import Chem
except ImportError:
    print("Error: RDKit not installed. Run: pip install rdkit")
    sys.exit(1)


def get_memory_mb():
    """Get current process memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def benchmark_sdf_parse(sdf_file: str, iterations: int = 5):
    """Benchmark SDF parsing with RDKit."""
    results = {
        "tool": "rdkit",
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

        # Use SDMolSupplier for streaming parsing
        supplier = Chem.SDMolSupplier(str(sdf_file))

        for mol in supplier:
            if mol is not None:
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


def benchmark_sdf_from_string(sdf_content: str, iterations: int = 100):
    """Benchmark parsing SDF from string (in-memory)."""
    results = {
        "tool": "rdkit",
        "type": "string_parse",
        "iterations": iterations,
    }

    # Warm-up
    for _ in range(10):
        Chem.MolFromMolBlock(sdf_content.split("$$$$")[0])

    start_time = time.perf_counter()
    for _ in range(iterations):
        # Parse single molecule from string
        mol_block = sdf_content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block)
        if mol is None:
            print("Warning: Failed to parse molecule")

    elapsed = time.perf_counter() - start_time
    results["elapsed_seconds"] = elapsed
    results["parses_per_second"] = iterations / elapsed
    results["ns_per_parse"] = (elapsed / iterations) * 1e9

    return results


def main():
    parser = argparse.ArgumentParser(description="Benchmark RDKit SDF parsing")
    parser.add_argument("sdf_file", help="Path to SDF file to benchmark")
    parser.add_argument("--iterations", "-n", type=int, default=5,
                        help="Number of benchmark iterations (default: 5)")
    parser.add_argument("--output", "-o", help="Output JSON file for results")
    args = parser.parse_args()

    sdf_path = Path(args.sdf_file)
    if not sdf_path.exists():
        print(f"Error: File not found: {sdf_path}")
        sys.exit(1)

    print(f"Benchmarking RDKit on: {sdf_path}")
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
