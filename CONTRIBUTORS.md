# Contributors

Thank you to everyone who has contributed to sdfrust!

## Maintainers

- **Hosein Fooladi** ([@HFooladi](https://github.com/HFooladi)) - Creator and lead maintainer

## Development

```bash
# Clone and build
git clone https://github.com/hfooladi/sdfrust.git
cd sdfrust
cargo build

# Run tests
cargo test                          # All Rust tests (~530)
cargo test --features geometry      # Include geometry-gated tests

# Linting and formatting
cargo clippy --workspace --features geometry
cargo fmt --check

# Benchmarks
cargo bench

# Python bindings
cd sdfrust-python
pip install maturin
maturin develop --features numpy
pytest tests/ -v                    # ~346 Python tests
```

### Project Structure

```
src/           Rust library (parsers, writers, descriptors, featurization)
tests/         Integration tests and test data
benches/       Criterion benchmarks
sdfrust-python/  Python bindings (PyO3) with its own tests and examples
```

See `CLAUDE.md` for a detailed module-level architecture overview.

## Contributing

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Adding New Functionality

- **Adding a feature:** Create module in `src/` → export from `lib.rs` → write tests in `tests/` → update `ROADMAP.md`
- **Adding a file format:** Create parser in `src/parser/<fmt>.rs` → map to `Molecule` → add integration tests with real files → add writer if needed

### Code Conventions

- Use `thiserror` for error types; return `Result<T, SdfError>` from fallible functions
- Accept `BufRead` trait for parser input (not concrete types)
- All public items require doc comments
- Internal indices are 0-based (convert to/from 1-based at parse/write boundaries)
- Unit and integration tests are expected for all new functionality

## Acknowledgments

- The cheminformatics community for SDF and MOL2 format specifications
- The Rust community for excellent tooling and libraries
