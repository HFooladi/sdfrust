//! Transparent gzip decompression support for file parsing.
//!
//! This module provides utilities for automatically handling gzip-compressed
//! files based on their `.gz` extension. It is only available when the `gzip`
//! feature is enabled.

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::error::Result;

/// A reader that transparently handles both plain and gzip-compressed files.
///
/// This enum wraps either a plain `BufReader<File>` or a `BufReader<GzDecoder<File>>`,
/// allowing the same parsing code to work with both compressed and uncompressed files.
pub enum MaybeGzReader {
    /// Plain (uncompressed) file reader.
    Plain(BufReader<File>),
    /// Gzip-compressed file reader.
    Gzip(BufReader<GzDecoder<File>>),
}

impl Read for MaybeGzReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            MaybeGzReader::Plain(reader) => reader.read(buf),
            MaybeGzReader::Gzip(reader) => reader.read(buf),
        }
    }
}

impl BufRead for MaybeGzReader {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        match self {
            MaybeGzReader::Plain(reader) => reader.fill_buf(),
            MaybeGzReader::Gzip(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            MaybeGzReader::Plain(reader) => reader.consume(amt),
            MaybeGzReader::Gzip(reader) => reader.consume(amt),
        }
    }
}

/// Checks if a path has a `.gz` extension.
///
/// # Example
///
/// ```ignore
/// use sdfrust::parser::compression::is_gzip_path;
///
/// assert!(is_gzip_path("molecule.sdf.gz"));
/// assert!(!is_gzip_path("molecule.sdf"));
/// ```
pub fn is_gzip_path<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

/// Opens a file and returns a reader that handles gzip decompression if needed.
///
/// If the file path ends in `.gz`, the file is decompressed transparently.
/// Otherwise, the file is read as-is.
///
/// # Example
///
/// ```ignore
/// use sdfrust::parser::compression::open_maybe_gz;
///
/// // Works with both compressed and uncompressed files
/// let reader = open_maybe_gz("molecule.sdf.gz")?;
/// let reader = open_maybe_gz("molecule.sdf")?;
/// ```
pub fn open_maybe_gz<P: AsRef<Path>>(path: P) -> Result<MaybeGzReader> {
    let file = File::open(&path)?;

    if is_gzip_path(&path) {
        let decoder = GzDecoder::new(file);
        Ok(MaybeGzReader::Gzip(BufReader::new(decoder)))
    } else {
        Ok(MaybeGzReader::Plain(BufReader::new(file)))
    }
}

/// Opens a gzip file and reads its entire contents into a String.
///
/// This is useful for format detection functions that need to read
/// the entire file into memory before parsing.
pub fn read_maybe_gz_to_string<P: AsRef<Path>>(path: P) -> Result<String> {
    let file = File::open(&path)?;
    let mut content = String::new();

    if is_gzip_path(&path) {
        let mut decoder = GzDecoder::new(file);
        decoder.read_to_string(&mut content)?;
    } else {
        let mut reader = BufReader::new(file);
        reader.read_to_string(&mut content)?;
    }

    Ok(content)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_gzip_path() {
        assert!(is_gzip_path("test.sdf.gz"));
        assert!(is_gzip_path("test.mol2.GZ"));
        assert!(is_gzip_path("test.xyz.Gz"));
        assert!(!is_gzip_path("test.sdf"));
        assert!(!is_gzip_path("test.mol2"));
        assert!(!is_gzip_path("test.gz.sdf")); // .gz not at end
    }
}
