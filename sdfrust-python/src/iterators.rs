//! Python bindings for streaming iterators over SDF and MOL2 files.

use pyo3::exceptions::PyStopIteration;
use pyo3::prelude::*;
use std::path::Path;

#[cfg(not(feature = "gzip"))]
use std::fs::File;
#[cfg(not(feature = "gzip"))]
use std::io::BufReader;

use sdfrust::{iter_mol2_file, iter_sdf_file, iter_sdf_v3000_file, iter_xyz_file};

#[cfg(not(feature = "gzip"))]
use sdfrust::{Mol2Iterator, SdfIterator, SdfV3000Iterator, XyzIterator};

#[cfg(feature = "gzip")]
use sdfrust::parser::MaybeGzReader;
#[cfg(feature = "gzip")]
use sdfrust::{Mol2Iterator, SdfIterator, SdfV3000Iterator, XyzIterator};

use crate::error::convert_error;
use crate::molecule::PyMolecule;

// ============================================================================
// SDF V2000 Iterator
// ============================================================================

/// Iterator over molecules in an SDF file (V2000 format).
///
/// This provides memory-efficient iteration over large SDF files
/// without loading all molecules into memory at once.
///
/// When compiled with the `gzip` feature, this iterator can also
/// read gzip-compressed files (`.sdf.gz`).
#[pyclass(name = "SdfIterator")]
pub struct PySdfIterator {
    #[cfg(not(feature = "gzip"))]
    inner: SdfIterator<BufReader<File>>,
    #[cfg(feature = "gzip")]
    inner: SdfIterator<MaybeGzReader>,
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
/// When compiled with the `gzip` feature, files ending in `.gz` are
/// automatically decompressed.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///     ValueError: If gzip file is provided but gzip feature is not enabled.
///
/// Example:
///     >>> for mol in iter_sdf_file("database.sdf"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
///     >>> # With gzip feature enabled:
///     >>> for mol in iter_sdf_file("database.sdf.gz"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "iter_sdf_file")]
pub fn py_iter_sdf_file(path: &str) -> PyResult<PySdfIterator> {
    let iter = iter_sdf_file(Path::new(path)).map_err(convert_error)?;
    Ok(PySdfIterator { inner: iter })
}

// ============================================================================
// SDF V3000 Iterator
// ============================================================================

/// Iterator over molecules in an SDF file (V3000 format).
///
/// When compiled with the `gzip` feature, this iterator can also
/// read gzip-compressed files (`.sdf.gz`).
#[pyclass(name = "SdfV3000Iterator")]
pub struct PySdfV3000Iterator {
    #[cfg(not(feature = "gzip"))]
    inner: SdfV3000Iterator<BufReader<File>>,
    #[cfg(feature = "gzip")]
    inner: SdfV3000Iterator<MaybeGzReader>,
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
/// When compiled with the `gzip` feature, files ending in `.gz` are
/// automatically decompressed.
///
/// Args:
///     path: Path to the SDF file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///     ValueError: If gzip file is provided but gzip feature is not enabled.
#[pyfunction]
#[pyo3(name = "iter_sdf_v3000_file")]
pub fn py_iter_sdf_v3000_file(path: &str) -> PyResult<PySdfV3000Iterator> {
    let iter = iter_sdf_v3000_file(Path::new(path)).map_err(convert_error)?;
    Ok(PySdfV3000Iterator { inner: iter })
}

// ============================================================================
// MOL2 Iterator
// ============================================================================

/// Iterator over molecules in a MOL2 file.
///
/// This provides memory-efficient iteration over large MOL2 files
/// without loading all molecules into memory at once.
///
/// When compiled with the `gzip` feature, this iterator can also
/// read gzip-compressed files (`.mol2.gz`).
#[pyclass(name = "Mol2Iterator")]
pub struct PyMol2Iterator {
    #[cfg(not(feature = "gzip"))]
    inner: Mol2Iterator<BufReader<File>>,
    #[cfg(feature = "gzip")]
    inner: Mol2Iterator<MaybeGzReader>,
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
/// When compiled with the `gzip` feature, files ending in `.gz` are
/// automatically decompressed.
///
/// Args:
///     path: Path to the MOL2 file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///     ValueError: If gzip file is provided but gzip feature is not enabled.
///
/// Example:
///     >>> for mol in iter_mol2_file("database.mol2"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
///     >>> # With gzip feature enabled:
///     >>> for mol in iter_mol2_file("database.mol2.gz"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "iter_mol2_file")]
pub fn py_iter_mol2_file(path: &str) -> PyResult<PyMol2Iterator> {
    let iter = iter_mol2_file(Path::new(path)).map_err(convert_error)?;
    Ok(PyMol2Iterator { inner: iter })
}

// ============================================================================
// XYZ Iterator
// ============================================================================

/// Iterator over molecules in an XYZ file.
///
/// This provides memory-efficient iteration over large XYZ files
/// without loading all molecules into memory at once.
///
/// When compiled with the `gzip` feature, this iterator can also
/// read gzip-compressed files (`.xyz.gz`).
#[pyclass(name = "XyzIterator")]
pub struct PyXyzIterator {
    #[cfg(not(feature = "gzip"))]
    inner: XyzIterator<BufReader<File>>,
    #[cfg(feature = "gzip")]
    inner: XyzIterator<MaybeGzReader>,
}

#[pymethods]
impl PyXyzIterator {
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

/// Create an iterator over molecules in an XYZ file.
///
/// This is memory-efficient for large files as molecules are parsed
/// one at a time rather than loading the entire file.
///
/// When compiled with the `gzip` feature, files ending in `.gz` are
/// automatically decompressed.
///
/// Args:
///     path: Path to the XYZ file.
///
/// Returns:
///     An iterator that yields Molecule objects.
///
/// Raises:
///     IOError: If the file cannot be opened.
///     ValueError: If gzip file is provided but gzip feature is not enabled.
///
/// Example:
///     >>> for mol in iter_xyz_file("trajectory.xyz"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
///     >>> # With gzip feature enabled:
///     >>> for mol in iter_xyz_file("trajectory.xyz.gz"):
///     ...     print(f"{mol.name}: {mol.num_atoms} atoms")
#[pyfunction]
#[pyo3(name = "iter_xyz_file")]
pub fn py_iter_xyz_file(path: &str) -> PyResult<PyXyzIterator> {
    let iter = iter_xyz_file(Path::new(path)).map_err(convert_error)?;
    Ok(PyXyzIterator { inner: iter })
}
