use thiserror::Error;

/// Errors that can occur when parsing or writing molecular structure files.
#[derive(Error, Debug)]
pub enum SdfError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error at line {line}: {message}")]
    Parse { line: usize, message: String },

    #[error("Invalid atom count: expected {expected}, found {found}")]
    AtomCountMismatch { expected: usize, found: usize },

    #[error("Invalid bond count: expected {expected}, found {found}")]
    BondCountMismatch { expected: usize, found: usize },

    #[error("Invalid atom index {index} in bond (molecule has {atom_count} atoms)")]
    InvalidAtomIndex { index: usize, atom_count: usize },

    #[error("Invalid bond order: {0}")]
    InvalidBondOrder(u8),

    #[error("Invalid counts line format: {0}")]
    InvalidCountsLine(String),

    #[error("Missing required section: {0}")]
    MissingSection(String),

    #[error("Empty file")]
    EmptyFile,

    #[error("Invalid coordinate value: {0}")]
    InvalidCoordinate(String),

    #[error("Invalid charge value: {0}")]
    InvalidCharge(String),
}

/// Result type alias for SDF operations.
pub type Result<T> = std::result::Result<T, SdfError>;
