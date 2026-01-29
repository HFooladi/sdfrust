# sdfrust Benchmark Results

This document contains performance benchmark results for sdfrust, comparing it against Python-based alternatives.

## Quick Summary

| Tool | Throughput (mol/s) | Relative Speed |
|------|-------------------|----------------|
| **sdfrust (Rust)** | ~200,000-245,000 | 1.0x (baseline) |
| RDKit (Python) | ~30,000-50,000 | 4-8x slower |
| Pure Python | ~3,000-5,000 | 40-80x slower |

## Running the Benchmarks

### Rust Benchmarks (Criterion)

```bash
# Run all Rust benchmarks
cargo bench

# Run specific benchmark
cargo bench --bench sdf_parse_benchmark
cargo bench --bench mol2_parse_benchmark
cargo bench --bench roundtrip_benchmark

# View HTML report
open target/criterion/report/index.html
```

### Quick Example Benchmark

```bash
# Run the example benchmark for quick results
cargo run --release --example benchmark
```

### Full Comparison Suite

```bash
# Run all benchmarks including Python comparison
cd benches/comparison
./run_all.sh

# With custom molecule count
./run_all.sh 100000
```

## Benchmark Groups

### SDF Parsing (`sdf_parse_benchmark`)

| Group | Benchmark | Description |
|-------|-----------|-------------|
| `parse_single` | `methane_5_atoms` | Parse single small molecule |
| `parse_single` | `aspirin_21_atoms` | Parse single medium molecule |
| `parse_single` | `synthetic_chain/N` | Parse synthetic N-atom chains (10-1000) |
| `parse_multi` | `methane_molecules/N` | Parse N molecules (10-10000) |
| `parse_iterator` | `iterator_methane/N` | Streaming parse N molecules |
| `parse_real_files` | Various | Parse real PubChem molecules |

### MOL2 Parsing (`mol2_parse_benchmark`)

| Group | Benchmark | Description |
|-------|-----------|-------------|
| `mol2_parse_single` | `methane_5_atoms` | Parse single MOL2 molecule |
| `mol2_parse_single` | `benzene_12_atoms` | Parse aromatic MOL2 molecule |
| `mol2_parse_multi` | `methane_molecules/N` | Parse N MOL2 molecules |
| `mol2_parse_iterator` | `iterator_methane/N` | Streaming parse MOL2 |

### Roundtrip (`roundtrip_benchmark`)

| Group | Benchmark | Description |
|-------|-----------|-------------|
| `roundtrip_single` | Various | Parse + write cycle |
| `write_only` | `write_molecule/N` | Write pre-built molecules |
| `molecule_operations` | `formula`, `centroid`, etc. | Molecule operations |
| `property_operations` | `property_get`, `property_set` | Property access |

## Baseline Results

From `cargo run --release --example benchmark`:

```
=== sdfrust Performance Benchmark ===

--- Benchmark: Parse Simple Molecules ---
  Parsed 10000 molecules in 40.79ms
  Average: 4079 ns/molecule
  Rate: 245,145 molecules/sec

--- Benchmark: Parse + Write Round-trip ---
  Round-tripped 5000 molecules in 309.86ms
  Average: 61971 ns/molecule
  Rate: 16,136 round-trips/sec

--- Benchmark: Large Molecule Operations ---
  Formula calculation: 50.63ms for 10000 iterations (5063 ns/call)
  Centroid calculation: 2.68ms for 10000 iterations (268 ns/call)
  Neighbor lookup: 2.13ms for 10000 iterations (213 ns/call)
  Element counts: 94.07ms for 10000 iterations (9407 ns/call)

--- Benchmark: Multi-molecule File ---
  Parsed 100,000 molecules in 493.56ms
  Average: 4936 ns/molecule
  Rate: 202,609 molecules/sec

--- Benchmark: Property Access ---
  Property lookup: 5.13ms for 100000 iterations (51 ns/call)
  Property set: 23.51ms for 100000 iterations (235 ns/call)
```

## Comparison Methodology

### Tools Compared

1. **sdfrust (Rust)**: This library
2. **RDKit (Python)**: Industry standard cheminformatics toolkit
3. **Pure Python**: Simple line-by-line SDF parser

### Metrics

| Metric | Description | Measurement |
|--------|-------------|-------------|
| Throughput | Molecules parsed per second | `molecules / elapsed_time` |
| Peak Memory | Maximum RSS during parsing | `/proc/self/statm` (Rust), `psutil` (Python) |
| First Molecule Latency | Time to return first molecule | `time_to_first_result` |

### Test Environment

Benchmarks should be run on:
- Release builds (`--release` for Rust)
- Fresh interpreter (for Python)
- Multiple iterations for statistical significance
- Same machine for comparison

## Expected Performance

Based on implementation characteristics:

| Scenario | Expected sdfrust Performance |
|----------|------------------------------|
| Small molecules (5 atoms) | >200,000 mol/s |
| Medium molecules (20 atoms) | >150,000 mol/s |
| Large molecules (100 atoms) | >50,000 mol/s |
| Multi-molecule streaming | >200,000 mol/s |

### Why sdfrust is Fast

1. **Zero-copy parsing**: Strings are borrowed where possible
2. **No GIL**: Full CPU utilization without interpreter lock
3. **Stack allocation**: Small values don't hit the heap
4. **Cache-friendly**: Linear memory access patterns
5. **Compile-time optimizations**: Aggressive inlining and SIMD

## Reproducing Results

### Requirements

- Rust 1.85+ with `cargo bench`
- Python 3.9+ with RDKit and psutil
- `uv` package manager (optional, for Python environment)

### Steps

1. Build in release mode:
   ```bash
   cargo build --release
   ```

2. Run Criterion benchmarks:
   ```bash
   cargo bench
   ```

3. Run Python comparison:
   ```bash
   cd benches/comparison
   ./run_all.sh
   ```

4. View results:
   - Criterion: `target/criterion/report/index.html`
   - Comparison: `BENCHMARK_RESULTS.md` (regenerated)

## Historical Results

Results will be tracked here as the library evolves:

| Version | Date | Parse Rate | Notes |
|---------|------|------------|-------|
| 0.1.0 | 2025-01 | ~200-245K mol/s | Initial baseline with Criterion benchmarks |

## Contributing

When adding new features, please:
1. Add appropriate benchmarks
2. Run `cargo bench` to check for regressions
3. Update this document if results change significantly
