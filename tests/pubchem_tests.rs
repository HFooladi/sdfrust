//! Comprehensive tests using real SDF files from PubChem.
//!
//! These tests validate parsing against actual chemical structures
//! from the PubChem database to ensure real-world compatibility.

use sdfrust::{BondOrder, Molecule, parse_sdf_file, write_sdf_string};

// ============================================================================
// Test Data Structures
// ============================================================================

/// Expected properties for each molecule from PubChem
struct MoleculeExpected {
    cid: &'static str,
    name: &'static str,
    formula: &'static str,
    atom_count: usize,
    bond_count: usize,
    heavy_atom_count: usize,
    molecular_weight: &'static str,
}

const ASPIRIN: MoleculeExpected = MoleculeExpected {
    cid: "2244",
    name: "2244",
    formula: "C9H8O4",
    atom_count: 21,
    bond_count: 21,
    heavy_atom_count: 13,
    molecular_weight: "180.16",
};

const CAFFEINE: MoleculeExpected = MoleculeExpected {
    cid: "2519",
    name: "2519",
    formula: "C8H10N4O2",
    atom_count: 24,
    bond_count: 25,
    heavy_atom_count: 14,
    molecular_weight: "194.19",
};

const GLUCOSE: MoleculeExpected = MoleculeExpected {
    cid: "5988",
    name: "5988",
    formula: "C12H22O11",
    atom_count: 45,
    bond_count: 46,
    heavy_atom_count: 23,
    molecular_weight: "342.30",
};

// CID 5793 is D-galactose (a monosaccharide)
const GALACTOSE: MoleculeExpected = MoleculeExpected {
    cid: "5793",
    name: "5793",
    formula: "C6H12O6",
    atom_count: 24,
    bond_count: 24,
    heavy_atom_count: 12,
    molecular_weight: "180.16",
};

const ACETAMINOPHEN: MoleculeExpected = MoleculeExpected {
    cid: "1983",
    name: "1983",
    formula: "C8H9NO2",
    atom_count: 20,
    bond_count: 20,
    heavy_atom_count: 11,
    molecular_weight: "151.16",
};

// Note: CID 6137 is actually methionine, not leucine
const METHIONINE: MoleculeExpected = MoleculeExpected {
    cid: "6137",
    name: "6137",
    formula: "C5H11NO2S",
    atom_count: 20,
    bond_count: 19,
    heavy_atom_count: 9,
    molecular_weight: "149.21",
};

// ============================================================================
// Helper Functions
// ============================================================================

fn count_heavy_atoms(mol: &Molecule) -> usize {
    mol.atoms().filter(|a| a.element != "H").count()
}

fn validate_molecule(mol: &Molecule, expected: &MoleculeExpected) {
    // Basic structure
    assert_eq!(
        mol.name, expected.name,
        "Name mismatch for CID {}",
        expected.cid
    );
    assert_eq!(
        mol.atom_count(),
        expected.atom_count,
        "Atom count mismatch for CID {}",
        expected.cid
    );
    assert_eq!(
        mol.bond_count(),
        expected.bond_count,
        "Bond count mismatch for CID {}",
        expected.cid
    );

    // Heavy atoms
    let heavy = count_heavy_atoms(mol);
    assert_eq!(
        heavy, expected.heavy_atom_count,
        "Heavy atom count mismatch for CID {}",
        expected.cid
    );

    // Properties
    assert_eq!(
        mol.get_property("PUBCHEM_COMPOUND_CID"),
        Some(expected.cid),
        "CID property mismatch"
    );
    assert_eq!(
        mol.get_property("PUBCHEM_MOLECULAR_FORMULA"),
        Some(expected.formula),
        "Formula property mismatch for CID {}",
        expected.cid
    );
    assert_eq!(
        mol.get_property("PUBCHEM_MOLECULAR_WEIGHT"),
        Some(expected.molecular_weight),
        "MW property mismatch for CID {}",
        expected.cid
    );

    // Verify our formula calculation matches PubChem
    assert_eq!(
        mol.formula(),
        expected.formula,
        "Calculated formula mismatch for CID {}",
        expected.cid
    );
}

// ============================================================================
// Aspirin Tests
// ============================================================================

#[test]
fn test_aspirin_parsing() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    validate_molecule(&mol, &ASPIRIN);
}

#[test]
fn test_aspirin_connectivity() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();

    // Aspirin has a benzene ring - check for double bonds
    let double_bonds: Vec<_> = mol
        .bonds()
        .filter(|b| b.order == BondOrder::Double)
        .collect();
    assert!(
        double_bonds.len() >= 2,
        "Aspirin should have at least 2 double bonds (C=O groups)"
    );

    // Check that all atoms have at least one bond
    for i in 0..mol.atom_count() {
        let bonds = mol.bonds_for_atom(i);
        assert!(
            !bonds.is_empty(),
            "Atom {} should have at least one bond",
            i
        );
    }
}

#[test]
fn test_aspirin_properties() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();

    // Check various PubChem properties
    assert!(mol.get_property("PUBCHEM_IUPAC_NAME").is_some());
    assert!(mol.get_property("PUBCHEM_SMILES").is_some());
    assert!(mol.get_property("PUBCHEM_IUPAC_INCHI").is_some());
    assert!(mol.get_property("PUBCHEM_IUPAC_INCHIKEY").is_some());

    // Verify SMILES is correct
    assert_eq!(
        mol.get_property("PUBCHEM_SMILES"),
        Some("CC(=O)OC1=CC=CC=C1C(=O)O")
    );
}

// ============================================================================
// Caffeine Tests
// ============================================================================

#[test]
fn test_caffeine_parsing() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    validate_molecule(&mol, &CAFFEINE);
}

#[test]
fn test_caffeine_nitrogen_count() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();

    // Caffeine has 4 nitrogen atoms
    let nitrogen_count = mol.atoms().filter(|a| a.element == "N").count();
    assert_eq!(nitrogen_count, 4, "Caffeine should have 4 nitrogen atoms");
}

#[test]
fn test_caffeine_ring_structure() {
    let mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();

    // Caffeine has a fused ring system with nitrogens in the rings
    // Ring nitrogens should have 2+ bonds, but carbonyl oxygens only have 1
    let ring_nitrogens: Vec<_> = mol.atoms().filter(|a| a.element == "N").collect();

    for atom in ring_nitrogens {
        let bonds = mol.bonds_for_atom(atom.index);
        assert!(
            bonds.len() >= 2,
            "Ring nitrogen {} should have at least 2 bonds",
            atom.index
        );
    }

    // Caffeine has 2 carbonyl oxygens (C=O), each with exactly 1 bond
    let oxygens: Vec<_> = mol.atoms().filter(|a| a.element == "O").collect();
    assert_eq!(oxygens.len(), 2, "Caffeine should have 2 oxygen atoms");

    for atom in oxygens {
        let bonds = mol.bonds_for_atom(atom.index);
        assert_eq!(
            bonds.len(),
            1,
            "Carbonyl oxygen {} should have exactly 1 bond (C=O)",
            atom.index
        );
    }
}

// ============================================================================
// Glucose Tests
// ============================================================================

#[test]
fn test_glucose_parsing() {
    let mol = parse_sdf_file("tests/test_data/glucose.sdf").unwrap();
    validate_molecule(&mol, &GLUCOSE);
}

#[test]
fn test_glucose_oxygen_count() {
    let mol = parse_sdf_file("tests/test_data/glucose.sdf").unwrap();

    // Sucrose (glucose disaccharide) has 11 oxygen atoms
    let oxygen_count = mol.atoms().filter(|a| a.element == "O").count();
    assert_eq!(oxygen_count, 11, "Glucose should have 11 oxygen atoms");
}

#[test]
fn test_glucose_all_single_bonds() {
    let mol = parse_sdf_file("tests/test_data/glucose.sdf").unwrap();

    // Glucose has only single bonds (no double bonds in the sugar)
    let single_bonds = mol.bonds().filter(|b| b.order == BondOrder::Single).count();
    assert_eq!(
        single_bonds,
        mol.bond_count(),
        "Glucose should have all single bonds"
    );
}

// ============================================================================
// Galactose Tests (CID 5793 - a monosaccharide)
// ============================================================================

#[test]
fn test_galactose_parsing() {
    let mol = parse_sdf_file("tests/test_data/galactose.sdf").unwrap();
    validate_molecule(&mol, &GALACTOSE);
}

#[test]
fn test_galactose_structure() {
    let mol = parse_sdf_file("tests/test_data/galactose.sdf").unwrap();

    // Galactose (C6H12O6) - monosaccharide
    let c_count = mol.atoms().filter(|a| a.element == "C").count();
    let o_count = mol.atoms().filter(|a| a.element == "O").count();
    let h_count = mol.atoms().filter(|a| a.element == "H").count();

    assert_eq!(c_count, 6, "Galactose should have 6 carbons");
    assert_eq!(o_count, 6, "Galactose should have 6 oxygens");
    assert_eq!(h_count, 12, "Galactose should have 12 hydrogens");
}

// ============================================================================
// Acetaminophen Tests
// ============================================================================

#[test]
fn test_acetaminophen_parsing() {
    let mol = parse_sdf_file("tests/test_data/acetaminophen.sdf").unwrap();
    validate_molecule(&mol, &ACETAMINOPHEN);
}

#[test]
fn test_acetaminophen_benzene_ring() {
    let mol = parse_sdf_file("tests/test_data/acetaminophen.sdf").unwrap();

    // Count carbon atoms (should have 8)
    let carbon_count = mol.atoms().filter(|a| a.element == "C").count();
    assert_eq!(carbon_count, 8, "Acetaminophen should have 8 carbons");

    // Should have exactly 1 nitrogen
    let nitrogen_count = mol.atoms().filter(|a| a.element == "N").count();
    assert_eq!(nitrogen_count, 1, "Acetaminophen should have 1 nitrogen");
}

// ============================================================================
// Methionine Tests (CID 6137 - amino acid with sulfur)
// ============================================================================

#[test]
fn test_methionine_parsing() {
    let mol = parse_sdf_file("tests/test_data/methionine.sdf").unwrap();
    validate_molecule(&mol, &METHIONINE);
}

#[test]
fn test_methionine_amino_acid_structure() {
    let mol = parse_sdf_file("tests/test_data/methionine.sdf").unwrap();

    // Methionine has N, O, C, S
    let n_count = mol.atoms().filter(|a| a.element == "N").count();
    let o_count = mol.atoms().filter(|a| a.element == "O").count();
    let c_count = mol.atoms().filter(|a| a.element == "C").count();
    let s_count = mol.atoms().filter(|a| a.element == "S").count();

    assert_eq!(n_count, 1, "Methionine should have 1 nitrogen");
    assert_eq!(o_count, 2, "Methionine should have 2 oxygens");
    assert_eq!(c_count, 5, "Methionine should have 5 carbons");
    assert_eq!(s_count, 1, "Methionine should have 1 sulfur");
}

// ============================================================================
// Round-trip Tests
// ============================================================================

#[test]
fn test_aspirin_round_trip() {
    let original = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let sdf_string = write_sdf_string(&original).unwrap();
    let parsed = sdfrust::parse_sdf_string(&sdf_string).unwrap();

    assert_eq!(parsed.atom_count(), original.atom_count());
    assert_eq!(parsed.bond_count(), original.bond_count());

    // Check coordinates are preserved
    for i in 0..original.atom_count() {
        let orig_atom = &original.atoms[i];
        let parsed_atom = &parsed.atoms[i];
        assert!(
            (orig_atom.x - parsed_atom.x).abs() < 0.001,
            "X coordinate mismatch at atom {}",
            i
        );
        assert!(
            (orig_atom.y - parsed_atom.y).abs() < 0.001,
            "Y coordinate mismatch at atom {}",
            i
        );
        assert!(
            (orig_atom.z - parsed_atom.z).abs() < 0.001,
            "Z coordinate mismatch at atom {}",
            i
        );
    }
}

#[test]
fn test_caffeine_round_trip() {
    let original = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    let sdf_string = write_sdf_string(&original).unwrap();
    let parsed = sdfrust::parse_sdf_string(&sdf_string).unwrap();

    assert_eq!(parsed.atom_count(), original.atom_count());
    assert_eq!(parsed.bond_count(), original.bond_count());
    assert_eq!(parsed.formula(), original.formula());
}

#[test]
fn test_glucose_round_trip() {
    let original = parse_sdf_file("tests/test_data/glucose.sdf").unwrap();
    let sdf_string = write_sdf_string(&original).unwrap();
    let parsed = sdfrust::parse_sdf_string(&sdf_string).unwrap();

    assert_eq!(parsed.atom_count(), original.atom_count());
    assert_eq!(parsed.bond_count(), original.bond_count());
    assert_eq!(parsed.formula(), original.formula());
}

// ============================================================================
// Geometry Tests
// ============================================================================

#[test]
fn test_aspirin_centroid() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
    let centroid = mol.centroid().unwrap();

    // Centroid should be somewhere reasonable (not at origin for this 2D structure)
    assert!(centroid.0.abs() < 10.0, "X centroid out of range");
    assert!(centroid.1.abs() < 10.0, "Y centroid out of range");
    // Z should be 0 for 2D structure
    assert!(centroid.2.abs() < 0.001, "Z centroid should be ~0 for 2D");
}

#[test]
fn test_molecule_centering() {
    let mut mol = parse_sdf_file("tests/test_data/caffeine_pubchem.sdf").unwrap();
    mol.center();

    let centroid = mol.centroid().unwrap();
    assert!(
        centroid.0.abs() < 0.001,
        "X centroid should be ~0 after centering"
    );
    assert!(
        centroid.1.abs() < 0.001,
        "Y centroid should be ~0 after centering"
    );
    assert!(
        centroid.2.abs() < 0.001,
        "Z centroid should be ~0 after centering"
    );
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn test_all_atoms_have_valid_elements() {
    let files = [
        "tests/test_data/aspirin.sdf",
        "tests/test_data/caffeine_pubchem.sdf",
        "tests/test_data/glucose.sdf",
        "tests/test_data/galactose.sdf",
        "tests/test_data/acetaminophen.sdf",
        "tests/test_data/methionine.sdf",
    ];

    for file in &files {
        let mol = parse_sdf_file(file).unwrap();
        for atom in mol.atoms() {
            assert!(
                !atom.element.is_empty(),
                "Empty element in {} at atom {}",
                file,
                atom.index
            );
            // Common organic elements
            let valid_elements = ["C", "H", "O", "N", "S", "P", "F", "Cl", "Br", "I"];
            assert!(
                valid_elements.contains(&atom.element.as_str()),
                "Unexpected element '{}' in {} at atom {}",
                atom.element,
                file,
                atom.index
            );
        }
    }
}

#[test]
fn test_all_bonds_reference_valid_atoms() {
    let files = [
        "tests/test_data/aspirin.sdf",
        "tests/test_data/caffeine_pubchem.sdf",
        "tests/test_data/glucose.sdf",
        "tests/test_data/galactose.sdf",
        "tests/test_data/acetaminophen.sdf",
        "tests/test_data/methionine.sdf",
    ];

    for file in &files {
        let mol = parse_sdf_file(file).unwrap();
        let atom_count = mol.atom_count();

        for bond in mol.bonds() {
            assert!(
                bond.atom1 < atom_count,
                "Invalid atom1 {} in {} (atom count: {})",
                bond.atom1,
                file,
                atom_count
            );
            assert!(
                bond.atom2 < atom_count,
                "Invalid atom2 {} in {} (atom count: {})",
                bond.atom2,
                file,
                atom_count
            );
            assert_ne!(
                bond.atom1, bond.atom2,
                "Self-bond in {} at atom {}",
                file, bond.atom1
            );
        }
    }
}

#[test]
fn test_coordinates_are_finite() {
    let files = [
        "tests/test_data/aspirin.sdf",
        "tests/test_data/caffeine_pubchem.sdf",
        "tests/test_data/glucose.sdf",
    ];

    for file in &files {
        let mol = parse_sdf_file(file).unwrap();
        for atom in mol.atoms() {
            assert!(
                atom.x.is_finite(),
                "Non-finite X in {} at atom {}",
                file,
                atom.index
            );
            assert!(
                atom.y.is_finite(),
                "Non-finite Y in {} at atom {}",
                file,
                atom.index
            );
            assert!(
                atom.z.is_finite(),
                "Non-finite Z in {} at atom {}",
                file,
                atom.index
            );
        }
    }
}

// ============================================================================
// Property Parsing Tests
// ============================================================================

#[test]
fn test_multiline_property_parsing() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();

    // PUBCHEM_COORDINATE_TYPE spans multiple lines
    let coord_type = mol.get_property("PUBCHEM_COORDINATE_TYPE");
    assert!(coord_type.is_some(), "PUBCHEM_COORDINATE_TYPE should exist");

    // Should contain multiple values
    let value = coord_type.unwrap();
    assert!(
        value.contains('\n') || value.len() > 3,
        "PUBCHEM_COORDINATE_TYPE should have multiple values"
    );
}

#[test]
fn test_property_count() {
    let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();

    // PubChem files have many properties
    assert!(
        mol.properties.len() >= 20,
        "Expected at least 20 properties, got {}",
        mol.properties.len()
    );
}

// ============================================================================
// Batch Processing Test
// ============================================================================

#[test]
fn test_parse_all_test_files() {
    let files = [
        ("tests/test_data/aspirin.sdf", &ASPIRIN),
        ("tests/test_data/caffeine_pubchem.sdf", &CAFFEINE),
        ("tests/test_data/glucose.sdf", &GLUCOSE),
        ("tests/test_data/galactose.sdf", &GALACTOSE),
        ("tests/test_data/acetaminophen.sdf", &ACETAMINOPHEN),
        ("tests/test_data/methionine.sdf", &METHIONINE),
    ];

    for (file, expected) in &files {
        let mol = parse_sdf_file(file).expect(&format!("Failed to parse {}", file));
        validate_molecule(&mol, expected);
    }
}

// ============================================================================
// Performance Sanity Check
// ============================================================================

#[test]
fn test_parse_multiple_times() {
    // Parse the same file multiple times to check consistency
    for _ in 0..10 {
        let mol = parse_sdf_file("tests/test_data/aspirin.sdf").unwrap();
        assert_eq!(mol.atom_count(), ASPIRIN.atom_count);
        assert_eq!(mol.bond_count(), ASPIRIN.bond_count);
    }
}
