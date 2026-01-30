//! Python bindings for streaming iterators over SDF and MOL2 files.

use pyo3::exceptions::PyStopIteration;
use pyo3::prelude::*;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use sdfrust::{
    iter_mol2_file, iter_sdf_file, iter_sdf_v3000_file, Mol2Iterator, SdfIterator, SdfV3000Iterator,
};

use crate::error::convert_error;
use crate::molecule::PyMolecule;

/// Iterator over molecules in an SDF file (V2000 format).
///
/// This provides memory-efficient iteration over large SDF files
/// without loading all molecules into memory at once.
#[pyclass(name = "SdfIterator")]
pub struct PySdfIterator {
    inner: SdfIterator<BufReader<File>>,
}

#[pymethods]
impl PySdfIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyMolecule>> {
        match slf.inner.next() {
            Some(Ok(mol)) => Ok(Some(PyMolecule::from(mol))),
            Some(Err(e)) => Err(convert_error(e)),
            None => Err(PyStopIteration::new_err(())),
        }
    }
}

/// Create an iterator over molecules in an SDF file (V2000 format).
///
/// This is memory-efficient for large files as molecules are parsed
/// one at a time rather than loading the entire file.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///
/// Example:
///     >>> for mol in iter_sdf_file("database.sdf"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "iter_sdf_file")]
pub fn py_iter_sdf_file(path: &str) -> PyResult<PySdfIterator> {
    let iter = iter_sdf_file(Path::new(path)).map_err(convert_error)?;
    Ok(PySdfIterator { inner: iter })
}

/// Iterator over molecules in an SDF file (V3000 format).
#[pyclass(name = "SdfV3000Iterator")]
pub struct PySdfV3000Iterator {
    inner: SdfV3000Iterator<BufReader<File>>,
}

#[pymethods]
impl PySdfV3000Iterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyMolecule>> {
        match slf.inner.next() {
            Some(Ok(mol)) => Ok(Some(PyMolecule::from(mol))),
            Some(Err(e)) => Err(convert_error(e)),
            None => Err(PyStopIteration::new_err(())),
        }
    }
}

/// Create an iterator over molecules in an SDF file (V3000 format).
///
/// This is memory-efficient for large files as molecules are parsed
/// one at a time rather than loading the entire file.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
#[pyfunction]
#[pyo3(name = "iter_sdf_v3000_file")]
pub fn py_iter_sdf_v3000_file(path: &str) -> PyResult<PySdfV3000Iterator> {
    let iter = iter_sdf_v3000_file(Path::new(path)).map_err(convert_error)?;
    Ok(PySdfV3000Iterator { inner: iter })
}

/// Iterator over molecules in a MOL2 file.
///
/// This provides memory-efficient iteration over large MOL2 files
/// without loading all molecules into memory at once.
#[pyclass(name = "Mol2Iterator")]
pub struct PyMol2Iterator {
    inner: Mol2Iterator<BufReader<File>>,
}

#[pymethods]
impl PyMol2Iterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> PyResult<Option<PyMolecule>> {
        match slf.inner.next() {
            Some(Ok(mol)) => Ok(Some(PyMolecule::from(mol))),
            Some(Err(e)) => Err(convert_error(e)),
            None => Err(PyStopIteration::new_err(())),
        }
    }
}

/// Create an iterator over molecules in a MOL2 file.
///
/// This is memory-efficient for large files as molecules are parsed
/// one at a time rather than loading the entire file.
///
/// Args:
///     path: Path to the MOL2 file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///
/// Example:
///     >>> for mol in iter_mol2_file("database.mol2"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "iter_mol2_file")]
pub fn py_iter_mol2_file(path: &str) -> PyResult<PyMol2Iterator> {
    let iter = iter_mol2_file(Path::new(path)).map_err(convert_error)?;
    Ok(PyMol2Iterator { inner: iter })
}
