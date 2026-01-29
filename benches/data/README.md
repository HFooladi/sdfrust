# Benchmark Test Data

This directory contains test data for sdfrust performance benchmarks.

## Synthetic Data Generation

Generate synthetic test SDF files using the provided script:

```bash
# Generate 10,000 molecules (default)
python generate_synthetic.py -n 10000 -o synthetic_10000.sdf

# Generate 100,000 molecules with varied sizes
python generate_synthetic.py -n 100000 --varied -o synthetic_100k_varied.sdf

# Generate small test file
python generate_synthetic.py -n 100 -o synthetic_100.sdf
```

## Real-World Datasets

For realistic benchmarks, you can download subsets from public chemical databases:

### PubChem Compound

Download from: https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/

```bash
# Download a single chunk (~25MB compressed, ~500K molecules)
wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/Compound_000000001_000500000.sdf.gz
gunzip Compound_000000001_000500000.sdf.gz
```

### ChEMBL

Download from: https://www.ebi.ac.uk/chembl/

```bash
# The full ChEMBL database (~2 million molecules)
wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_33.sdf.gz
gunzip chembl_33.sdf.gz
```

### ZINC Database

Download from: https://zinc.docking.org/

ZINC provides subsets for drug-like compounds in various formats.

## Recommended Test Files

| Test Type | File | Molecules | Purpose |
|-----------|------|-----------|---------|
| Quick test | synthetic_100.sdf | 100 | Fast iteration during development |
| Standard | synthetic_10000.sdf | 10,000 | Default benchmark size |
| Large | synthetic_100000.sdf | 100,000 | Stress test / realistic load |
| Real-world | PubChem chunk | ~500,000 | Real molecule diversity |

## File Sizes (Approximate)

- 1,000 molecules: ~100 KB
- 10,000 molecules: ~1 MB
- 100,000 molecules: ~10 MB
- 1,000,000 molecules: ~100 MB

## Usage in Benchmarks

The benchmark scripts automatically look for test data in this directory:

```bash
# Run all benchmarks with default synthetic data
../comparison/run_all.sh

# Run with specific molecule count
../comparison/run_all.sh 100000
```

## Gitignore

Generated SDF files are ignored by git to avoid bloating the repository.
Only keep small test files (<1MB) in version control.
