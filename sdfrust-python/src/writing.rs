//! Python bindings for SDF writing functions.

use pyo3::prelude::*;
use std::path::Path;

use sdfrust::{
    write_sdf_auto_file, write_sdf_auto_string, write_sdf_file, write_sdf_file_multi,
    write_sdf_string, write_sdf_v3000_file, write_sdf_v3000_file_multi, write_sdf_v3000_string,
};

use crate::error::convert_error;
use crate::molecule::PyMolecule;

// ============================================================
// SDF V2000 Writing
// ============================================================

/// Write a molecule to an SDF file (V2000 format).
///
/// Args:
///     molecule: The molecule to write.
///     path: Path to the output file.
///
/// Raises:
///     IOError: If the file cannot be written.
#[pyfunction]
#[pyo3(name = "write_sdf_file")]
pub fn py_write_sdf_file(molecule: &PyMolecule, path: &str) -> PyResult<()> {
    write_sdf_file(Path::new(path), &molecule.inner).map_err(convert_error)
}

/// Write a molecule to an SDF string (V2000 format).
///
/// Args:
///     molecule: The molecule to write.
///
/// Returns:
///     The SDF content as a string.
#[pyfunction]
#[pyo3(name = "write_sdf_string")]
pub fn py_write_sdf_string(molecule: &PyMolecule) -> PyResult<String> {
    write_sdf_string(&molecule.inner).map_err(convert_error)
}

/// Write multiple molecules to an SDF file (V2000 format).
///
/// Args:
///     molecules: List of molecules to write.
///     path: Path to the output file.
///
/// Raises:
///     IOError: If the file cannot be written.
#[pyfunction]
#[pyo3(name = "write_sdf_file_multi")]
pub fn py_write_sdf_file_multi(molecules: Vec<PyMolecule>, path: &str) -> PyResult<()> {
    let mols: Vec<_> = molecules.into_iter().map(|m| m.inner).collect();
    write_sdf_file_multi(Path::new(path), &mols).map_err(convert_error)
}

// ============================================================
// SDF V3000 Writing
// ============================================================

/// Write a molecule to an SDF file (V3000 format).
///
/// Args:
///     molecule: The molecule to write.
///     path: Path to the output file.
///
/// Raises:
///     IOError: If the file cannot be written.
#[pyfunction]
#[pyo3(name = "write_sdf_v3000_file")]
pub fn py_write_sdf_v3000_file(molecule: &PyMolecule, path: &str) -> PyResult<()> {
    write_sdf_v3000_file(Path::new(path), &molecule.inner).map_err(convert_error)
}

/// Write a molecule to an SDF string (V3000 format).
///
/// Args:
///     molecule: The molecule to write.
///
/// Returns:
///     The SDF content as a string.
#[pyfunction]
#[pyo3(name = "write_sdf_v3000_string")]
pub fn py_write_sdf_v3000_string(molecule: &PyMolecule) -> PyResult<String> {
    write_sdf_v3000_string(&molecule.inner).map_err(convert_error)
}

/// Write multiple molecules to an SDF file (V3000 format).
///
/// Args:
///     molecules: List of molecules to write.
///     path: Path to the output file.
///
/// Raises:
///     IOError: If the file cannot be written.
#[pyfunction]
#[pyo3(name = "write_sdf_v3000_file_multi")]
pub fn py_write_sdf_v3000_file_multi(molecules: Vec<PyMolecule>, path: &str) -> PyResult<()> {
    let mols: Vec<_> = molecules.into_iter().map(|m| m.inner).collect();
    write_sdf_v3000_file_multi(Path::new(path), &mols).map_err(convert_error)
}

// ============================================================
// SDF Auto-Format Writing (V2000/V3000)
// ============================================================

/// Write a molecule to an SDF file with automatic format selection.
///
/// Uses V3000 format if the molecule has >999 atoms/bonds or V3000-only features,
/// otherwise uses V2000 format.
///
/// Args:
///     molecule: The molecule to write.
///     path: Path to the output file.
///
/// Raises:
///     IOError: If the file cannot be written.
#[pyfunction]
#[pyo3(name = "write_sdf_auto_file")]
pub fn py_write_sdf_auto_file(molecule: &PyMolecule, path: &str) -> PyResult<()> {
    write_sdf_auto_file(Path::new(path), &molecule.inner).map_err(convert_error)
}

/// Write a molecule to an SDF string with automatic format selection.
///
/// Uses V3000 format if the molecule has >999 atoms/bonds or V3000-only features,
/// otherwise uses V2000 format.
///
/// Args:
///     molecule: The molecule to write.
///
/// Returns:
///     The SDF content as a string.
#[pyfunction]
#[pyo3(name = "write_sdf_auto_string")]
pub fn py_write_sdf_auto_string(molecule: &PyMolecule) -> PyResult<String> {
    write_sdf_auto_string(&molecule.inner).map_err(convert_error)
}
