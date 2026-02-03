//! Python bindings for SDF and MOL2 parsing functions.

use pyo3::prelude::*;
use std::path::Path;

use sdfrust::{
    detect_format, parse_auto_file, parse_auto_file_multi, parse_auto_string,
    parse_auto_string_multi, parse_mol2_file, parse_mol2_file_multi, parse_mol2_string,
    parse_mol2_string_multi, parse_sdf_auto_file, parse_sdf_auto_file_multi, parse_sdf_auto_string,
    parse_sdf_auto_string_multi, parse_sdf_file, parse_sdf_file_multi, parse_sdf_string,
    parse_sdf_string_multi, parse_sdf_v3000_file, parse_sdf_v3000_file_multi,
    parse_sdf_v3000_string, parse_sdf_v3000_string_multi, parse_xyz_file, parse_xyz_file_multi,
    parse_xyz_string, parse_xyz_string_multi,
};

use crate::error::convert_error;
use crate::molecule::PyMolecule;

// ============================================================
// SDF V2000 Parsing
// ============================================================

/// Parse a single molecule from an SDF file (V2000 format).
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_file")]
pub fn py_parse_sdf_file(path: &str) -> PyResult<PyMolecule> {
    parse_sdf_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from an SDF string (V2000 format).
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_string")]
pub fn py_parse_sdf_string(content: &str) -> PyResult<PyMolecule> {
    parse_sdf_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF file (V2000 format).
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_file_multi")]
pub fn py_parse_sdf_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF string (V2000 format).
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_string_multi")]
pub fn py_parse_sdf_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

// ============================================================
// SDF V3000 Parsing
// ============================================================

/// Parse a single molecule from an SDF file (V3000 format).
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_v3000_file")]
pub fn py_parse_sdf_v3000_file(path: &str) -> PyResult<PyMolecule> {
    parse_sdf_v3000_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from an SDF string (V3000 format).
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_v3000_string")]
pub fn py_parse_sdf_v3000_string(content: &str) -> PyResult<PyMolecule> {
    parse_sdf_v3000_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF file (V3000 format).
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_v3000_file_multi")]
pub fn py_parse_sdf_v3000_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_v3000_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF string (V3000 format).
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_v3000_string_multi")]
pub fn py_parse_sdf_v3000_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_v3000_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

// ============================================================
// SDF Auto-Detection Parsing (V2000/V3000)
// ============================================================

/// Parse a single molecule from an SDF file with automatic format detection.
///
/// Automatically detects whether the file is V2000 or V3000 format.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_auto_file")]
pub fn py_parse_sdf_auto_file(path: &str) -> PyResult<PyMolecule> {
    parse_sdf_auto_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from an SDF string with automatic format detection.
///
/// Automatically detects whether the content is V2000 or V3000 format.
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_auto_string")]
pub fn py_parse_sdf_auto_string(content: &str) -> PyResult<PyMolecule> {
    parse_sdf_auto_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF file with automatic format detection.
///
/// Automatically detects whether the file is V2000 or V3000 format.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_auto_file_multi")]
pub fn py_parse_sdf_auto_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_auto_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from an SDF string with automatic format detection.
///
/// Automatically detects whether the content is V2000 or V3000 format.
///
/// Args:
///     content: The SDF content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_sdf_auto_string_multi")]
pub fn py_parse_sdf_auto_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_sdf_auto_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

// ============================================================
// MOL2 Parsing
// ============================================================

/// Parse a single molecule from a MOL2 file.
///
/// Args:
///     path: Path to the MOL2 file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_mol2_file")]
pub fn py_parse_mol2_file(path: &str) -> PyResult<PyMolecule> {
    parse_mol2_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from a MOL2 string.
///
/// Args:
///     content: The MOL2 content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_mol2_string")]
pub fn py_parse_mol2_string(content: &str) -> PyResult<PyMolecule> {
    parse_mol2_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from a MOL2 file.
///
/// Args:
///     path: Path to the MOL2 file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_mol2_file_multi")]
pub fn py_parse_mol2_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_mol2_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from a MOL2 string.
///
/// Args:
///     content: The MOL2 content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_mol2_string_multi")]
pub fn py_parse_mol2_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_mol2_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

// ============================================================
// Unified Auto-Detection (SDF V2000, V3000, MOL2)
// ============================================================

/// Detect the format of molecular structure file content.
///
/// This function examines the content to determine whether it is:
/// - MOL2 format (returns "mol2")
/// - SDF V3000 format (returns "sdf_v3000")
/// - SDF V2000 format (returns "sdf_v2000")
///
/// Args:
///     content: The file content as a string.
///
/// Returns:
///     A string indicating the detected format: "sdf_v2000", "sdf_v3000", or "mol2".
///
/// Example:
///     >>> import sdfrust
///     >>> format = sdfrust.detect_format("@<TRIPOS>MOLECULE\\ntest\\n")
///     >>> print(format)  # "mol2"
#[pyfunction]
#[pyo3(name = "detect_format")]
pub fn py_detect_format(content: &str) -> String {
    let format = detect_format(content);
    format.to_string()
}

/// Parse a single molecule from a file with automatic format detection.
///
/// This function reads the file, detects whether it is SDF V2000, V3000,
/// or MOL2 format, and uses the appropriate parser.
///
/// Args:
///     path: Path to the molecular structure file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
///
/// Example:
///     >>> import sdfrust
///     >>> mol = sdfrust.parse_auto_file("molecule.sdf")  # Works with .sdf or .mol2
///     >>> print(mol.name)
#[pyfunction]
#[pyo3(name = "parse_auto_file")]
pub fn py_parse_auto_file(path: &str) -> PyResult<PyMolecule> {
    parse_auto_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from a string with automatic format detection.
///
/// This function automatically detects whether the content is SDF V2000, V3000,
/// or MOL2 format and uses the appropriate parser.
///
/// Args:
///     content: The file content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_auto_string")]
pub fn py_parse_auto_string(content: &str) -> PyResult<PyMolecule> {
    parse_auto_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from a file with automatic format detection.
///
/// This function reads the file, detects whether it is SDF V2000, V3000,
/// or MOL2 format, and uses the appropriate parser.
///
/// Args:
///     path: Path to the molecular structure file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_auto_file_multi")]
pub fn py_parse_auto_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_auto_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from a string with automatic format detection.
///
/// This function automatically detects whether the content is SDF V2000, V3000,
/// or MOL2 format and uses the appropriate parser.
///
/// Args:
///     content: The file content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_auto_string_multi")]
pub fn py_parse_auto_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_auto_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

// ============================================================
// XYZ Parsing
// ============================================================

/// Parse a single molecule from an XYZ file.
///
/// XYZ is a simple format containing only atomic coordinates (no bonds).
///
/// Args:
///     path: Path to the XYZ file.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
///
/// Example:
///     >>> import sdfrust
///     >>> mol = sdfrust.parse_xyz_file("water.xyz")
///     >>> print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "parse_xyz_file")]
pub fn py_parse_xyz_file(path: &str) -> PyResult<PyMolecule> {
    parse_xyz_file(Path::new(path))
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse a single molecule from an XYZ string.
///
/// Args:
///     content: The XYZ content as a string.
///
/// Returns:
///     The parsed Molecule.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_xyz_string")]
pub fn py_parse_xyz_string(content: &str) -> PyResult<PyMolecule> {
    parse_xyz_string(content)
        .map(PyMolecule::from)
        .map_err(convert_error)
}

/// Parse multiple molecules from an XYZ file.
///
/// XYZ files can contain multiple molecules concatenated together.
///
/// Args:
///     path: Path to the XYZ file.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     IOError: If the file cannot be read.
///     ValueError: If the file cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_xyz_file_multi")]
pub fn py_parse_xyz_file_multi(path: &str) -> PyResult<Vec<PyMolecule>> {
    parse_xyz_file_multi(Path::new(path))
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}

/// Parse multiple molecules from an XYZ string.
///
/// Args:
///     content: The XYZ content as a string.
///
/// Returns:
///     A list of parsed Molecules.
///
/// Raises:
///     ValueError: If the content cannot be parsed.
#[pyfunction]
#[pyo3(name = "parse_xyz_string_multi")]
pub fn py_parse_xyz_string_multi(content: &str) -> PyResult<Vec<PyMolecule>> {
    parse_xyz_string_multi(content)
        .map(|mols| mols.into_iter().map(PyMolecule::from).collect())
        .map_err(convert_error)
}
