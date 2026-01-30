"""
sdfrust - Fast Rust-based SDF and MOL2 molecular structure file parser.

This package provides high-performance parsing and writing of SDF (Structure Data File)
and MOL2 (TRIPOS) molecular structure formats. It is implemented in Rust for speed
and safety.

Quick Start:
    >>> import sdfrust
    >>> mol = sdfrust.parse_sdf_file("molecule.sdf")
    >>> print(f"Name: {mol.name}")
    >>> print(f"Atoms: {mol.num_atoms}")
    >>> print(f"Formula: {mol.formula()}")

Features:
    - Parse SDF V2000 and V3000 formats
    - Parse TRIPOS MOL2 format
    - Write SDF V2000 and V3000 formats
    - Memory-efficient streaming iterators for large files
    - Molecular descriptors (MW, exact mass, ring count, etc.)
    - NumPy integration for coordinate arrays (optional)
"""

from sdfrust._sdfrust import (
    # Version
    __version__,
    # Classes
    Atom,
    Bond,
    BondOrder,
    BondStereo,
    Molecule,
    SdfFormat,
    SdfIterator,
    SdfV3000Iterator,
    Mol2Iterator,
    # SDF V2000 parsing
    parse_sdf_file,
    parse_sdf_string,
    parse_sdf_file_multi,
    parse_sdf_string_multi,
    # SDF V3000 parsing
    parse_sdf_v3000_file,
    parse_sdf_v3000_string,
    parse_sdf_v3000_file_multi,
    parse_sdf_v3000_string_multi,
    # SDF auto-detection parsing
    parse_sdf_auto_file,
    parse_sdf_auto_string,
    parse_sdf_auto_file_multi,
    parse_sdf_auto_string_multi,
    # MOL2 parsing
    parse_mol2_file,
    parse_mol2_string,
    parse_mol2_file_multi,
    parse_mol2_string_multi,
    # SDF V2000 writing
    write_sdf_file,
    write_sdf_string,
    write_sdf_file_multi,
    # SDF V3000 writing
    write_sdf_v3000_file,
    write_sdf_v3000_string,
    write_sdf_v3000_file_multi,
    # SDF auto-format writing
    write_sdf_auto_file,
    write_sdf_auto_string,
    # Iterators
    iter_sdf_file,
    iter_sdf_v3000_file,
    iter_mol2_file,
)

__all__ = [
    # Version
    "__version__",
    # Classes
    "Atom",
    "Bond",
    "BondOrder",
    "BondStereo",
    "Molecule",
    "SdfFormat",
    "SdfIterator",
    "SdfV3000Iterator",
    "Mol2Iterator",
    # SDF V2000 parsing
    "parse_sdf_file",
    "parse_sdf_string",
    "parse_sdf_file_multi",
    "parse_sdf_string_multi",
    # SDF V3000 parsing
    "parse_sdf_v3000_file",
    "parse_sdf_v3000_string",
    "parse_sdf_v3000_file_multi",
    "parse_sdf_v3000_string_multi",
    # SDF auto-detection parsing
    "parse_sdf_auto_file",
    "parse_sdf_auto_string",
    "parse_sdf_auto_file_multi",
    "parse_sdf_auto_string_multi",
    # MOL2 parsing
    "parse_mol2_file",
    "parse_mol2_string",
    "parse_mol2_file_multi",
    "parse_mol2_string_multi",
    # SDF V2000 writing
    "write_sdf_file",
    "write_sdf_string",
    "write_sdf_file_multi",
    # SDF V3000 writing
    "write_sdf_v3000_file",
    "write_sdf_v3000_string",
    "write_sdf_v3000_file_multi",
    # SDF auto-format writing
    "write_sdf_auto_file",
    "write_sdf_auto_string",
    # Iterators
    "iter_sdf_file",
    "iter_sdf_v3000_file",
    "iter_mol2_file",
]
