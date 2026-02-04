//! Integration tests for gzip decompression support.
//!
//! These tests verify that gzipped files are correctly decompressed
//! when parsed using the various file parsing functions.

#![cfg(feature = "gzip")]

use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::Write;
use tempfile::NamedTempFile;

// Test data
const SIMPLE_SDF: &str = r#"methane
  test    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"#;

const SIMPLE_MOL2: &str = r#"@<TRIPOS>MOLECULE
methane
 5 4 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 H1          0.6289    0.6289    0.6289 H         1 MOL       0.0000
      3 H2         -0.6289   -0.6289    0.6289 H         1 MOL       0.0000
      4 H3         -0.6289    0.6289   -0.6289 H         1 MOL       0.0000
      5 H4          0.6289   -0.6289   -0.6289 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
"#;

const SIMPLE_XYZ: &str = r#"3
water molecule
O  0.000000  0.000000  0.117300
H  0.756950  0.000000 -0.469200
H -0.756950  0.000000 -0.469200
"#;

const V3000_SDF: &str = r#"methane
  sdfrust   01012500003D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 H 0.6289 0.6289 0.6289 0
M  V30 3 H -0.6289 -0.6289 0.6289 0
M  V30 4 H -0.6289 0.6289 -0.6289 0
M  V30 5 H 0.6289 -0.6289 -0.6289 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
"#;

/// Helper to create a gzipped temp file with the given content
fn create_gzip_file(content: &str, extension: &str) -> NamedTempFile {
    let mut file = tempfile::Builder::new()
        .suffix(extension)
        .tempfile()
        .unwrap();

    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(content.as_bytes()).unwrap();
    let compressed = encoder.finish().unwrap();

    file.write_all(&compressed).unwrap();
    file.flush().unwrap();
    file
}

// ============================================================================
// SDF V2000 Tests
// ============================================================================

#[test]
fn test_parse_gzipped_sdf_file() {
    let file = create_gzip_file(SIMPLE_SDF, ".sdf.gz");
    let mol = sdfrust::parse_sdf_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
    assert_eq!(mol.formula(), "CH4");
}

#[test]
fn test_parse_gzipped_sdf_file_multi() {
    let multi_sdf = format!("{}{}", SIMPLE_SDF, SIMPLE_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let mols = sdfrust::parse_sdf_file_multi(file.path()).unwrap();

    assert_eq!(mols.len(), 2);
    assert_eq!(mols[0].name, "methane");
    assert_eq!(mols[1].name, "methane");
}

#[test]
fn test_iter_gzipped_sdf_file() {
    let multi_sdf = format!("{}{}{}", SIMPLE_SDF, SIMPLE_SDF, SIMPLE_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let iter = sdfrust::iter_sdf_file(file.path()).unwrap();

    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
    assert_eq!(mols.len(), 3);
}

// ============================================================================
// SDF V3000 Tests
// ============================================================================

#[test]
fn test_parse_gzipped_sdf_v3000_file() {
    let file = create_gzip_file(V3000_SDF, ".sdf.gz");
    let mol = sdfrust::parse_sdf_v3000_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_parse_gzipped_sdf_v3000_file_multi() {
    let multi_sdf = format!("{}{}", V3000_SDF, V3000_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let mols = sdfrust::parse_sdf_v3000_file_multi(file.path()).unwrap();

    assert_eq!(mols.len(), 2);
}

#[test]
fn test_iter_gzipped_sdf_v3000_file() {
    let multi_sdf = format!("{}{}", V3000_SDF, V3000_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let iter = sdfrust::iter_sdf_v3000_file(file.path()).unwrap();

    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
    assert_eq!(mols.len(), 2);
}

// ============================================================================
// MOL2 Tests
// ============================================================================

#[test]
fn test_parse_gzipped_mol2_file() {
    let file = create_gzip_file(SIMPLE_MOL2, ".mol2.gz");
    let mol = sdfrust::parse_mol2_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
    assert_eq!(mol.bond_count(), 4);
}

#[test]
fn test_parse_gzipped_mol2_file_multi() {
    let multi_mol2 = format!("{}{}", SIMPLE_MOL2, SIMPLE_MOL2);
    let file = create_gzip_file(&multi_mol2, ".mol2.gz");
    let mols = sdfrust::parse_mol2_file_multi(file.path()).unwrap();

    assert_eq!(mols.len(), 2);
}

#[test]
fn test_iter_gzipped_mol2_file() {
    let multi_mol2 = format!("{}{}{}", SIMPLE_MOL2, SIMPLE_MOL2, SIMPLE_MOL2);
    let file = create_gzip_file(&multi_mol2, ".mol2.gz");
    let iter = sdfrust::iter_mol2_file(file.path()).unwrap();

    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
    assert_eq!(mols.len(), 3);
}

// ============================================================================
// XYZ Tests
// ============================================================================

#[test]
fn test_parse_gzipped_xyz_file() {
    let file = create_gzip_file(SIMPLE_XYZ, ".xyz.gz");
    let mol = sdfrust::parse_xyz_file(file.path()).unwrap();

    assert_eq!(mol.name, "water molecule");
    assert_eq!(mol.atom_count(), 3);
    assert_eq!(mol.bond_count(), 0); // XYZ has no bonds
    assert_eq!(mol.formula(), "H2O");
}

#[test]
fn test_parse_gzipped_xyz_file_multi() {
    let multi_xyz = format!("{}{}", SIMPLE_XYZ, SIMPLE_XYZ);
    let file = create_gzip_file(&multi_xyz, ".xyz.gz");
    let mols = sdfrust::parse_xyz_file_multi(file.path()).unwrap();

    assert_eq!(mols.len(), 2);
}

#[test]
fn test_iter_gzipped_xyz_file() {
    let multi_xyz = format!("{}{}{}", SIMPLE_XYZ, SIMPLE_XYZ, SIMPLE_XYZ);
    let file = create_gzip_file(&multi_xyz, ".xyz.gz");
    let iter = sdfrust::iter_xyz_file(file.path()).unwrap();

    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
    assert_eq!(mols.len(), 3);
}

// ============================================================================
// Auto-detection Tests
// ============================================================================

#[test]
fn test_parse_auto_gzipped_sdf() {
    let file = create_gzip_file(SIMPLE_SDF, ".sdf.gz");
    let mol = sdfrust::parse_auto_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.formula(), "CH4");
}

#[test]
fn test_parse_auto_gzipped_v3000() {
    let file = create_gzip_file(V3000_SDF, ".sdf.gz");
    let mol = sdfrust::parse_auto_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
    assert_eq!(mol.atom_count(), 5);
}

#[test]
fn test_parse_auto_gzipped_mol2() {
    let file = create_gzip_file(SIMPLE_MOL2, ".mol2.gz");
    let mol = sdfrust::parse_auto_file(file.path()).unwrap();

    assert_eq!(mol.name, "methane");
}

#[test]
fn test_parse_auto_gzipped_xyz() {
    let file = create_gzip_file(SIMPLE_XYZ, ".xyz.gz");
    let mol = sdfrust::parse_auto_file(file.path()).unwrap();

    assert_eq!(mol.name, "water molecule");
    assert_eq!(mol.formula(), "H2O");
}

#[test]
fn test_parse_auto_file_multi_gzipped() {
    let multi_sdf = format!("{}{}", SIMPLE_SDF, SIMPLE_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let mols = sdfrust::parse_auto_file_multi(file.path()).unwrap();

    assert_eq!(mols.len(), 2);
}

#[test]
fn test_iter_auto_gzipped() {
    let multi_sdf = format!("{}{}{}", SIMPLE_SDF, SIMPLE_SDF, SIMPLE_SDF);
    let file = create_gzip_file(&multi_sdf, ".sdf.gz");
    let iter = sdfrust::iter_auto_file(file.path()).unwrap();

    let mols: Vec<_> = iter.map(|r| r.unwrap()).collect();
    assert_eq!(mols.len(), 3);
}

// ============================================================================
// Compression module Tests
// ============================================================================

#[test]
fn test_is_gzip_path() {
    use sdfrust::parser::is_gzip_path;

    assert!(is_gzip_path("test.sdf.gz"));
    assert!(is_gzip_path("test.mol2.GZ"));
    assert!(is_gzip_path("test.xyz.Gz"));
    assert!(!is_gzip_path("test.sdf"));
    assert!(!is_gzip_path("test.mol2"));
    assert!(!is_gzip_path("test.gz.sdf")); // .gz not at end
}

#[test]
fn test_open_maybe_gz_plain_file() {
    use sdfrust::parser::open_maybe_gz;
    use std::io::BufRead;

    // Create a plain (non-gzipped) temp file
    let mut file = tempfile::Builder::new().suffix(".sdf").tempfile().unwrap();
    file.write_all(SIMPLE_SDF.as_bytes()).unwrap();
    file.flush().unwrap();

    let mut reader = open_maybe_gz(file.path()).unwrap();
    let mut first_line = String::new();
    reader.read_line(&mut first_line).unwrap();

    assert_eq!(first_line.trim(), "methane");
}

#[test]
fn test_open_maybe_gz_gzipped_file() {
    use sdfrust::parser::open_maybe_gz;
    use std::io::BufRead;

    let file = create_gzip_file(SIMPLE_SDF, ".sdf.gz");

    let mut reader = open_maybe_gz(file.path()).unwrap();
    let mut first_line = String::new();
    reader.read_line(&mut first_line).unwrap();

    assert_eq!(first_line.trim(), "methane");
}

// ============================================================================
// Edge Case Tests
// ============================================================================

#[test]
fn test_gzip_case_insensitive_extension() {
    // Test various case combinations of .gz extension
    for ext in &[".sdf.gz", ".sdf.GZ", ".sdf.Gz", ".sdf.gZ"] {
        let file = create_gzip_file(SIMPLE_SDF, ext);
        let mol = sdfrust::parse_sdf_file(file.path()).unwrap();
        assert_eq!(mol.name, "methane");
    }
}

#[test]
fn test_empty_gzipped_file() {
    let file = create_gzip_file("", ".sdf.gz");
    let result = sdfrust::parse_sdf_file(file.path());

    assert!(result.is_err());
}

#[test]
fn test_plain_file_still_works_with_gzip_feature() {
    // Create a plain (non-gzipped) temp file to ensure we didn't break
    // normal file parsing when the gzip feature is enabled
    let mut file = tempfile::Builder::new().suffix(".sdf").tempfile().unwrap();
    file.write_all(SIMPLE_SDF.as_bytes()).unwrap();
    file.flush().unwrap();

    let mol = sdfrust::parse_sdf_file(file.path()).unwrap();
    assert_eq!(mol.name, "methane");
    assert_eq!(mol.formula(), "CH4");
}
