# sdfrust

A fast, pure-Rust parser for SDF (Structure Data File) and MOL2 chemical structure files.

## Features

- Parse SDF V2000 format files (single and multi-molecule)
- Write SDF V2000 format files
- Iterate over large SDF files without loading everything into memory
- Access atom coordinates, bonds, and molecule properties
- Zero external dependencies for parsing (only `thiserror` for error handling)

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
sdfrust = "0.1"
```

## Quick Start

### Parse a single molecule

```rust
use sdfrust::parse_sdf_file;

let molecule = parse_sdf_file("molecule.sdf")?;
println!("Name: {}", molecule.name);
println!("Atoms: {}", molecule.atom_count());
println!("Formula: {}", molecule.formula());
```

### Parse multiple molecules

```rust
use sdfrust::parse_sdf_file_multi;

let molecules = parse_sdf_file_multi("database.sdf")?;
for mol in &molecules {
    println!("{}: {} atoms", mol.name, mol.atom_count());
}
```

### Iterate over a large file (memory efficient)

```rust
use sdfrust::iter_sdf_file;

for result in iter_sdf_file("large_database.sdf")? {
    let mol = result?;
    // Process each molecule without loading all into memory
}
```

### Parse from string

```rust
use sdfrust::parse_sdf_string;

let sdf_content = "...";  // SDF file content
let mol = parse_sdf_string(sdf_content)?;
```

### Write molecules

```rust
use sdfrust::{Molecule, Atom, Bond, BondOrder, write_sdf_file};

let mut mol = Molecule::new("water");
mol.atoms.push(Atom::new(0, "O", 0.0, 0.0, 0.0));
mol.atoms.push(Atom::new(1, "H", 0.96, 0.0, 0.0));
mol.atoms.push(Atom::new(2, "H", -0.24, 0.93, 0.0));
mol.bonds.push(Bond::new(0, 1, BondOrder::Single));
mol.bonds.push(Bond::new(0, 2, BondOrder::Single));

write_sdf_file("water.sdf", &mol)?;
```

## Data Model

### Molecule

The main container for molecular data:

```rust
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub properties: HashMap<String, String>,
}
```

### Atom

Represents an atom with 3D coordinates:

```rust
pub struct Atom {
    pub index: usize,
    pub element: String,
    pub x: f64, pub y: f64, pub z: f64,
    pub formal_charge: i8,
    pub mass_difference: i8,
    pub stereo_parity: Option<u8>,
}
```

### Bond

Represents a bond between two atoms:

```rust
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub order: BondOrder,
    pub stereo: BondStereo,
}
```

### BondOrder

```rust
pub enum BondOrder {
    Single, Double, Triple, Aromatic,
    SingleOrDouble, SingleOrAromatic, DoubleOrAromatic, Any,
}
```

## Molecule Methods

```rust
// Atom/bond counts
mol.atom_count()
mol.bond_count()
mol.is_empty()

// Iterators
mol.atoms()
mol.bonds()

// Connectivity
mol.neighbors(atom_index)      // Get connected atom indices
mol.bonds_for_atom(atom_index) // Get bonds for an atom

// Properties
mol.get_property("MW")
mol.set_property("SMILES", "CCO")

// Chemistry
mol.formula()         // "C2H6O"
mol.total_charge()    // Sum of formal charges
mol.element_counts()  // HashMap of element counts
mol.has_aromatic_bonds()
mol.has_charges()

// Geometry
mol.centroid()        // Geometric center
mol.center()          // Move centroid to origin
mol.translate(x, y, z)

// Filtering
mol.atoms_by_element("C")
mol.bonds_by_order(BondOrder::Double)
```

## Supported Format

Currently supports **SDF V2000** format, which is the most widely used format for chemical structure data.

### Planned

- SDF V3000 format
- MOL2 (TRIPOS) format

## Performance

The library is designed for high performance:

- Streaming parser for memory-efficient processing of large files
- Minimal allocations during parsing
- Zero-copy where possible

## License

MIT License
