//! Integration tests for XYZ format parsing.

use sdfrust::{
    FileFormat, SdfError, detect_format, iter_auto_file, iter_xyz_file, parse_auto_file,
    parse_auto_string, parse_xyz_file, parse_xyz_file_multi, parse_xyz_string,
    parse_xyz_string_multi,
};

// ============================================================================
// Basic Parsing Tests
// ============================================================================

#[test]
fn test_parse_water_xyz_file() {
    let mol = parse_xyz_file("tests/test_data/water.xyz").unwrap();

    assert_eq!(mol.name, "water molecule");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 0); // XYZ has no bonds
    assert_eq!(mol.formula(), "H2O");

    // Check oxygen coordinates
    let oxygen = &mol.atoms[0];
    assert_eq!(oxygen.element, "O");
    assert!((oxygen.x - 0.0).abs() < 1e-6);
    assert!((oxygen.y - 0.0).abs() < 1e-6);
    assert!((oxygen.z - 0.1173).abs() < 1e-6);

    // Check hydrogen atoms
    assert_eq!(mol.atoms[1].element, "H");
    assert_eq!(mol.atoms[2].element, "H");
}

#[test]
fn test_parse_xyz_string() {
    let xyz = r#"5
methane
C   0.000000   0.000000   0.000000
H   0.628900   0.628900   0.628900
H  -0.628900  -0.628900   0.628900
H  -0.628900   0.628900  -0.628900
H   0.628900  -0.628900  -0.628900
"#;

    let mol = parse_xyz_string(xyz).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.formula(), "CH4");

    // Check carbon is at origin
    let carbon = &mol.atoms[0];
    assert_eq!(carbon.element, "C");
    assert!((carbon.x - 0.0).abs() < 1e-6);
    assert!((carbon.y - 0.0).abs() < 1e-6);
    assert!((carbon.z - 0.0).abs() < 1e-6);
}

// ============================================================================
// Multi-molecule Tests
// ============================================================================

#[test]
fn test_parse_multi_xyz_file() {
    let mols = parse_xyz_file_multi("tests/test_data/multi.xyz").unwrap();

    assert_eq!(mols.len(), 3);

    // First molecule: water
    assert_eq!(mols[0].name, "water molecule");
    assert_eq!(mols[0].atom_count(), 3);
    assert_eq!(mols[0].formula(), "H2O");

    // Second molecule: methane
    assert_eq!(mols[1].name, "methane");
    assert_eq!(mols[1].atom_count(), 5);
    assert_eq!(mols[1].formula(), "CH4");

    // Third molecule: H2
    assert_eq!(mols[2].name, "diatomic hydrogen");
    assert_eq!(mols[2].atom_count(), 2);
    assert_eq!(mols[2].formula(), "H2");
}

#[test]
fn test_parse_xyz_string_multi() {
    let xyz = r#"2
molecule 1
C  0.0  0.0  0.0
O  1.2  0.0  0.0
3
molecule 2
N  0.0  0.0  0.0
H  1.0  0.0  0.0
H -1.0  0.0  0.0
"#;

    let mols = parse_xyz_string_multi(xyz).unwrap();

    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "molecule 1");
    assert_eq!(mols[0].formula(), "CO");
    assert_eq!(mols[1].name, "molecule 2");
    assert_eq!(mols[1].formula(), "H2N");
}

// ============================================================================
// Iterator Tests
// ============================================================================

#[test]
fn test_xyz_iterator() {
    let iter = iter_xyz_file("tests/test_data/multi.xyz").unwrap();
    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();

    assert_eq!(mols.len(), 3);
    assert_eq!(mols[0].name, "water molecule");
    assert_eq!(mols[1].name, "methane");
    assert_eq!(mols[2].name, "diatomic hydrogen");
}

#[test]
fn test_xyz_iterator_count() {
    let iter = iter_xyz_file("tests/test_data/multi.xyz").unwrap();
    assert_eq!(iter.count(), 3);
}

// ============================================================================
// Auto-detection Tests
// ============================================================================

#[test]
fn test_detect_xyz_format() {
    let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n";
    assert_eq!(detect_format(xyz), FileFormat::Xyz);
}

#[test]
fn test_detect_xyz_with_atomic_numbers() {
    let xyz = "3\ntest\n8 0.0 0.0 0.0\n1 1.0 0.0 0.0\n1 -1.0 0.0 0.0\n";
    assert_eq!(detect_format(xyz), FileFormat::Xyz);
}

#[test]
fn test_parse_auto_string_xyz() {
    let xyz = r#"3
water
O  0.000000  0.000000  0.117300
H  0.756950  0.000000 -0.469200
H -0.756950  0.000000 -0.469200
"#;

    let mol = parse_auto_string(xyz).unwrap();
    assert_eq!(mol.name, "water");
    assert_eq!(mol.atom_count(), 3);
}

#[test]
fn test_parse_auto_file_xyz() {
    let mol = parse_auto_file("tests/test_data/water.xyz").unwrap();
    assert_eq!(mol.name, "water molecule");
    assert_eq!(mol.formula(), "H2O");
}

#[test]
fn test_iter_auto_file_xyz() {
    let iter = iter_auto_file("tests/test_data/multi.xyz").unwrap();
    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();

    assert_eq!(mols.len(), 3);
}

// ============================================================================
// Edge Cases
// ============================================================================

#[test]
fn test_xyz_with_atomic_numbers() {
    let xyz = r#"3
atomic number test
8  0.0  0.0  0.0
1  1.0  0.0  0.0
1 -1.0  0.0  0.0
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert_eq!(mol.atoms[0].element, "O"); // 8 -> O
    assert_eq!(mol.atoms[1].element, "H"); // 1 -> H
    assert_eq!(mol.atoms[2].element, "H"); // 1 -> H
}

#[test]
fn test_xyz_case_normalization() {
    let xyz = r#"3
case test
o  0.0  0.0  0.0
h  1.0  0.0  0.0
CA -1.0  0.0  0.0
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert_eq!(mol.atoms[0].element, "O"); // o -> O
    assert_eq!(mol.atoms[1].element, "H"); // h -> H
    assert_eq!(mol.atoms[2].element, "Ca"); // CA -> Ca
}

#[test]
fn test_xyz_extra_columns_ignored() {
    let xyz = r#"2
extra columns
C  0.0  0.0  0.0  0.5  charge  extra_info
H  1.0  0.0  0.0 -0.2
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.atoms[0].element, "C");
    assert_eq!(mol.atoms[1].element, "H");
}

#[test]
fn test_xyz_with_blank_lines_between_molecules() {
    let xyz = r#"2
mol1
C  0.0  0.0  0.0
H  1.0  0.0  0.0

2
mol2
N  0.0  0.0  0.0
O  1.0  0.0  0.0
"#;

    let mols = parse_xyz_string_multi(xyz).unwrap();
    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "mol1");
    assert_eq!(mols[1].name, "mol2");
}

#[test]
fn test_xyz_negative_coordinates() {
    let xyz = r#"1
negative coords
C  -1.234567  -8.901234  -5.678901
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert!((mol.atoms[0].x - (-1.234567)).abs() < 1e-6);
    assert!((mol.atoms[0].y - (-8.901234)).abs() < 1e-6);
    assert!((mol.atoms[0].z - (-5.678901)).abs() < 1e-6);
}

#[test]
fn test_xyz_scientific_notation() {
    let xyz = r#"1
scientific notation
C  1.0e-5  2.5E+3  -3.7e2
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert!((mol.atoms[0].x - 1.0e-5).abs() < 1e-10);
    assert!((mol.atoms[0].y - 2.5e3).abs() < 1e-6);
    assert!((mol.atoms[0].z - (-3.7e2)).abs() < 1e-6);
}

// ============================================================================
// Error Cases
// ============================================================================

#[test]
fn test_xyz_empty_file() {
    let result = parse_xyz_string("");
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), SdfError::EmptyFile));
}

#[test]
fn test_xyz_only_whitespace() {
    let result = parse_xyz_string("   \n\n  \n");
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), SdfError::EmptyFile));
}

#[test]
fn test_xyz_invalid_atom_count() {
    let xyz = "abc\ntest\nC 0.0 0.0 0.0\n";
    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        SdfError::InvalidCountsLine(_)
    ));
}

#[test]
fn test_xyz_negative_atom_count() {
    // Negative numbers can't be parsed as usize
    let xyz = "-1\ntest\n";
    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
}

#[test]
fn test_xyz_fewer_atoms_than_declared() {
    let xyz = r#"5
missing atoms
C  0.0  0.0  0.0
H  1.0  0.0  0.0
"#;

    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
    match result.unwrap_err() {
        SdfError::AtomCountMismatch { expected, found } => {
            assert_eq!(expected, 5);
            assert_eq!(found, 2);
        }
        _ => panic!("Expected AtomCountMismatch error"),
    }
}

#[test]
fn test_xyz_invalid_coordinate() {
    let xyz = r#"1
bad coords
C  abc  0.0  0.0
"#;

    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
    assert!(matches!(
        result.unwrap_err(),
        SdfError::InvalidCoordinate(_)
    ));
}

#[test]
fn test_xyz_too_few_columns() {
    let xyz = r#"1
too few
C  0.0  0.0
"#;

    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), SdfError::Parse { .. }));
}

#[test]
fn test_xyz_missing_comment_line() {
    // Only atom count, no comment line
    let xyz = "1\n";
    let result = parse_xyz_string(xyz);
    assert!(result.is_err());
}

#[test]
fn test_xyz_file_not_found() {
    let result = parse_xyz_file("nonexistent_file.xyz");
    assert!(result.is_err());
    assert!(matches!(result.unwrap_err(), SdfError::Io(_)));
}

// ============================================================================
// Format Detection Edge Cases
// ============================================================================

#[test]
fn test_detect_sdf_not_confused_with_xyz() {
    // SDF content should not be detected as XYZ
    let sdf = r#"methane
  test    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
"#;

    assert_eq!(detect_format(sdf), FileFormat::SdfV2000);
}

#[test]
fn test_detect_mol2_not_confused_with_xyz() {
    let mol2 = "@<TRIPOS>MOLECULE\ntest\n";
    assert_eq!(detect_format(mol2), FileFormat::Mol2);
}

// ============================================================================
// Molecule Properties
// ============================================================================

#[test]
fn test_xyz_molecule_has_no_bonds() {
    let mol = parse_xyz_file("tests/test_data/water.xyz").unwrap();
    assert!(mol.bonds.is_empty());
    assert!(!mol.has_aromatic_bonds());
}

#[test]
fn test_xyz_molecule_centroid() {
    let xyz = r#"2
diatomic
C  0.0  0.0  0.0
C  2.0  0.0  0.0
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    let (cx, cy, cz) = mol.centroid().unwrap();
    assert!((cx - 1.0).abs() < 1e-6);
    assert!((cy - 0.0).abs() < 1e-6);
    assert!((cz - 0.0).abs() < 1e-6);
}

#[test]
fn test_xyz_heavy_atom_count() {
    let mol = parse_xyz_file("tests/test_data/water.xyz").unwrap();
    assert_eq!(mol.heavy_atom_count(), 1); // Only oxygen
}

#[test]
fn test_xyz_molecular_formula_ordering() {
    // Carbon should come first, then hydrogen, then alphabetical
    let xyz = r#"5
formula test
O  0.0  0.0  0.0
C  1.0  0.0  0.0
N  2.0  0.0  0.0
H  3.0  0.0  0.0
H  4.0  0.0  0.0
"#;

    let mol = parse_xyz_string(xyz).unwrap();
    assert_eq!(mol.formula(), "CH2NO"); // C first, H second, then N, O alphabetically
}
