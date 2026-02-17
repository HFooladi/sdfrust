//! Python bindings for sdfrust - a fast SDF, MOL2, and XYZ molecular structure parser.
//!
//! This crate provides Python bindings for the sdfrust library using PyO3.
//! It exposes all the core functionality for parsing and writing molecular
//! structure files in SDF, MOL2, and XYZ formats.

use pyo3::prelude::*;

pub mod atom;
pub mod bond;
pub mod error;
pub mod iterators;
pub mod molecule;
pub mod parsing;
pub mod writing;

/// Check if gzip support is enabled.
///
/// Returns True if the library was compiled with gzip support,
/// allowing parsing of `.gz` compressed files.
///
/// Example:
///     >>> import sdfrust
///     >>> if sdfrust.gzip_enabled():
///     ...     mol = sdfrust.parse_sdf_file("molecule.sdf.gz")
///     ... else:
///     ...     print("Gzip support not available")
#[pyfunction]
#[pyo3(name = "gzip_enabled")]
fn py_gzip_enabled() -> bool {
    cfg!(feature = "gzip")
}

use atom::PyAtom;
use bond::{PyBond, PyBondOrder, PyBondStereo};
use iterators::{PyMol2Iterator, PySdfIterator, PySdfV3000Iterator, PyXyzIterator};
use molecule::{PyMolecule, PySdfFormat};

/// sdfrust - Fast Rust-based SDF, MOL2, and XYZ molecular structure file parser.
///
/// This module provides high-performance parsing and writing of SDF (Structure Data File)
/// and MOL2 (TRIPOS) molecular structure formats. It is implemented in Rust for speed
/// and safety.
///
/// Quick Start:
///     >>> import sdfrust
///     >>> mol = sdfrust.parse_sdf_file("molecule.sdf")
///     >>> print(f"Name: {mol.name}")
///     >>> print(f"Atoms: {mol.num_atoms}")
///     >>> print(f"Formula: {mol.formula()}")
///
/// Features:
///     - Parse SDF V2000 and V3000 formats
///     - Parse TRIPOS MOL2 format
///     - Parse XYZ format
///     - Write SDF V2000 and V3000 formats
///     - Memory-efficient streaming iterators for large files
///     - Molecular descriptors (MW, exact mass, ring count, etc.)
///     - NumPy integration for coordinate arrays (optional)
///     - Transparent gzip decompression (optional, requires gzip feature)
#[pymodule]
fn _sdfrust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Version
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    // Classes
    m.add_class::<PyAtom>()?;
    m.add_class::<PyBond>()?;
    m.add_class::<PyBondOrder>()?;
    m.add_class::<PyBondStereo>()?;
    m.add_class::<PyMolecule>()?;
    m.add_class::<PySdfFormat>()?;
    m.add_class::<PySdfIterator>()?;
    m.add_class::<PySdfV3000Iterator>()?;
    m.add_class::<PyMol2Iterator>()?;
    m.add_class::<PyXyzIterator>()?;

    // SDF V2000 parsing functions
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_string_multi, m)?)?;

    // SDF V3000 parsing functions
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_v3000_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_v3000_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_v3000_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(
        parsing::py_parse_sdf_v3000_string_multi,
        m
    )?)?;

    // SDF auto-detection parsing functions (V2000/V3000 only)
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_auto_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_auto_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_sdf_auto_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(
        parsing::py_parse_sdf_auto_string_multi,
        m
    )?)?;

    // Unified auto-detection functions (SDF V2000, V3000, and MOL2)
    m.add_function(wrap_pyfunction!(parsing::py_detect_format, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_auto_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_auto_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_auto_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_auto_string_multi, m)?)?;

    // MOL2 parsing functions
    m.add_function(wrap_pyfunction!(parsing::py_parse_mol2_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_mol2_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_mol2_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_mol2_string_multi, m)?)?;

    // XYZ parsing functions
    m.add_function(wrap_pyfunction!(parsing::py_parse_xyz_file, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_xyz_string, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_xyz_file_multi, m)?)?;
    m.add_function(wrap_pyfunction!(parsing::py_parse_xyz_string_multi, m)?)?;

    // SDF V2000 writing functions
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_file, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_string, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_file_multi, m)?)?;

    // SDF V3000 writing functions
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_v3000_file, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_v3000_string, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_v3000_file_multi, m)?)?;

    // SDF auto-format writing functions
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_auto_file, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_sdf_auto_string, m)?)?;

    // MOL2 writing functions
    m.add_function(wrap_pyfunction!(writing::py_write_mol2_file, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_mol2_string, m)?)?;
    m.add_function(wrap_pyfunction!(writing::py_write_mol2_file_multi, m)?)?;

    // Iterator functions
    m.add_function(wrap_pyfunction!(iterators::py_iter_sdf_file, m)?)?;
    m.add_function(wrap_pyfunction!(iterators::py_iter_sdf_v3000_file, m)?)?;
    m.add_function(wrap_pyfunction!(iterators::py_iter_mol2_file, m)?)?;
    m.add_function(wrap_pyfunction!(iterators::py_iter_xyz_file, m)?)?;

    // Utility functions
    m.add_function(wrap_pyfunction!(py_gzip_enabled, m)?)?;

    Ok(())
}
