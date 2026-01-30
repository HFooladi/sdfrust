//! Python bindings for SDF and MOL2 parsing functions.

use pyo3::prelude::*;
use std::path::Path;

use sdfrust::{
    parse_mol2_file, parse_mol2_file_multi, parse_mol2_string, parse_mol2_string_multi,
    parse_sdf_auto_file, parse_sdf_auto_file_multi, parse_sdf_auto_string,
    parse_sdf_auto_string_multi, parse_sdf_file, parse_sdf_file_multi, parse_sdf_string,
    parse_sdf_string_multi, parse_sdf_v3000_file, parse_sdf_v3000_file_multi,
    parse_sdf_v3000_string, parse_sdf_v3000_string_multi,
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
