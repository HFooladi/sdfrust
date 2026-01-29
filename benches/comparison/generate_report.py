#!/usr/bin/env python3
"""
Generate a markdown comparison report from benchmark results.

Usage:
    python generate_report.py --rust <rust_results.json> --rdkit <rdkit_results.json> --python <python_results.json> --output BENCHMARK_RESULTS.md
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


def load_results(filepath: str) -> dict:
    """Load benchmark results from JSON file."""
    with open(filepath) as f:
        return json.load(f)


def format_number(n: float, decimals: int = 2) -> str:
    """Format a number with thousands separators."""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.{decimals}f}M"
    elif n >= 1_000:
        return f"{n / 1_000:.{decimals}f}K"
    else:
        return f"{n:.{decimals}f}"


def generate_report(
    rust_results: dict,
    rdkit_results: dict,
    python_results: dict,
    output_file: str,
):
    """Generate markdown comparison report."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    # Extract key metrics
    rust_throughput = rust_results.get("average", {}).get("molecules_per_second", 0)
    rdkit_throughput = rdkit_results.get("average", {}).get("molecules_per_second", 0)
    python_throughput = python_results.get("average", {}).get("molecules_per_second", 0)

    rust_memory = rust_results.get("average", {}).get("peak_memory_mb", 0)
    rdkit_memory = rdkit_results.get("average", {}).get("peak_memory_mb", 0)
    python_memory = python_results.get("average", {}).get("peak_memory_mb", 0)

    rust_latency = rust_results.get("average", {}).get("first_molecule_latency_ms", 0)
    rdkit_latency = rdkit_results.get("average", {}).get("first_molecule_latency_ms", 0)
    python_latency = python_results.get("average", {}).get("first_molecule_latency_ms", 0)

    # Calculate speedups
    rdkit_speedup = rust_throughput / rdkit_throughput if rdkit_throughput > 0 else 0
    python_speedup = rust_throughput / python_throughput if python_throughput > 0 else 0

    rdkit_memory_ratio = rdkit_memory / rust_memory if rust_memory > 0 else 0
    python_memory_ratio = python_memory / rust_memory if rust_memory > 0 else 0

    report = f"""# sdfrust Benchmark Results

Generated: {now}

## Summary

| Tool | Throughput (mol/s) | Peak Memory (MB) | First Mol Latency (ms) |
|------|-------------------|------------------|------------------------|
| **sdfrust (Rust)** | {format_number(rust_throughput)} | {rust_memory:.1f} | {rust_latency:.2f} |
| RDKit (Python) | {format_number(rdkit_throughput)} | {rdkit_memory:.1f} | {rdkit_latency:.2f} |
| Pure Python | {format_number(python_throughput)} | {python_memory:.1f} | {python_latency:.2f} |

## Performance Comparison

### Throughput

- **sdfrust is {rdkit_speedup:.1f}x faster than RDKit**
- **sdfrust is {python_speedup:.1f}x faster than pure Python**

### Memory Usage

- sdfrust uses {rdkit_memory_ratio:.1f}x less memory than RDKit
- sdfrust uses {python_memory_ratio:.1f}x less memory than pure Python

## Test Configuration

- **Test file**: {rust_results.get('file', 'N/A')}
- **Molecules**: {rust_results.get('runs', [{}])[0].get('molecules', 'N/A')}
- **Iterations**: {rust_results.get('iterations', 'N/A')}

## Detailed Results

### sdfrust (Rust)

```json
{json.dumps(rust_results.get('average', {}), indent=2)}
```

### RDKit (Python bindings)

```json
{json.dumps(rdkit_results.get('average', {}), indent=2)}
```

### Pure Python

```json
{json.dumps(python_results.get('average', {}), indent=2)}
```

## Methodology

All benchmarks were run on the same machine with:
- Warm-up iterations before measurement
- Multiple runs averaged
- Memory measured using /proc/self/statm (Rust) and psutil (Python)
- First molecule latency measures time to parse and return first molecule

### Rust Benchmarks

Run with Criterion:
```bash
cargo bench
```

### Python Benchmarks

Run with:
```bash
python benchmark_rdkit.py <sdf_file>
python benchmark_pure_python.py <sdf_file>
```

## Conclusions

sdfrust provides significant performance advantages over Python-based alternatives:

1. **Throughput**: {rdkit_speedup:.0f}-{python_speedup:.0f}x faster parsing
2. **Memory**: Lower memory footprint enables processing larger datasets
3. **Latency**: Fast first-molecule access for interactive use cases

These results validate sdfrust as a high-performance choice for SDF file processing
in performance-critical applications.
"""

    with open(output_file, "w") as f:
        f.write(report)

    print(f"Report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark comparison report")
    parser.add_argument("--rust", required=True, help="Rust benchmark results JSON")
    parser.add_argument("--rdkit", required=True, help="RDKit benchmark results JSON")
    parser.add_argument("--python", required=True, help="Pure Python benchmark results JSON")
    parser.add_argument("--output", "-o", default="BENCHMARK_RESULTS.md",
                        help="Output markdown file (default: BENCHMARK_RESULTS.md)")
    args = parser.parse_args()

    # Load results
    try:
        rust_results = load_results(args.rust)
    except FileNotFoundError:
        print(f"Error: Rust results not found: {args.rust}")
        sys.exit(1)

    try:
        rdkit_results = load_results(args.rdkit)
    except FileNotFoundError:
        print(f"Warning: RDKit results not found: {args.rdkit}")
        rdkit_results = {"average": {}}

    try:
        python_results = load_results(args.python)
    except FileNotFoundError:
        print(f"Warning: Python results not found: {args.python}")
        python_results = {"average": {}}

    generate_report(rust_results, rdkit_results, python_results, args.output)


if __name__ == "__main__":
    main()
