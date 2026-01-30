//! Error conversion from Rust SdfError to Python exceptions.

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use sdfrust::SdfError;

/// Convert SdfError to appropriate Python exception.
pub fn convert_error(err: SdfError) -> PyErr {
    match err {
        SdfError::Io(io_err) => PyIOError::new_err(io_err.to_string()),
        SdfError::Parse { line, message } => {
            PyValueError::new_err(format!("Parse error at line {}: {}", line, message))
        }
        SdfError::AtomCountMismatch { expected, found } => PyValueError::new_err(format!(
            "Invalid atom count: expected {}, found {}",
            expected, found
        )),
        SdfError::BondCountMismatch { expected, found } => PyValueError::new_err(format!(
            "Invalid bond count: expected {}, found {}",
            expected, found
        )),
        SdfError::InvalidAtomIndex { index, atom_count } => PyValueError::new_err(format!(
            "Invalid atom index {} in bond (molecule has {} atoms)",
            index, atom_count
        )),
        SdfError::InvalidBondOrder(order) => {
            PyValueError::new_err(format!("Invalid bond order: {}", order))
        }
        SdfError::InvalidCountsLine(line) => {
            PyValueError::new_err(format!("Invalid counts line format: {}", line))
        }
        SdfError::MissingSection(section) => {
            PyValueError::new_err(format!("Missing required section: {}", section))
        }
        SdfError::EmptyFile => PyValueError::new_err("Empty file"),
        SdfError::InvalidCoordinate(coord) => {
            PyValueError::new_err(format!("Invalid coordinate value: {}", coord))
        }
        SdfError::InvalidCharge(charge) => {
            PyValueError::new_err(format!("Invalid charge value: {}", charge))
        }
        SdfError::InvalidV3000Block(block) => {
            PyValueError::new_err(format!("Invalid V3000 block: {}", block))
        }
        SdfError::InvalidV3000AtomLine { line, message } => PyValueError::new_err(format!(
            "Invalid V3000 atom line at line {}: {}",
            line, message
        )),
        SdfError::InvalidV3000BondLine { line, message } => PyValueError::new_err(format!(
            "Invalid V3000 bond line at line {}: {}",
            line, message
        )),
        SdfError::AtomIdNotFound { id } => {
            PyValueError::new_err(format!("Atom ID {} not found in V3000 ID mapping", id))
        }
        SdfError::UnsupportedV3000Feature(feature) => {
            PyValueError::new_err(format!("Unsupported V3000 feature: {}", feature))
        }
    }
}
