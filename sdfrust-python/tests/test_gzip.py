"""Tests for gzip decompression support in sdfrust Python bindings."""

import gzip
import os
import tempfile

import pytest

import sdfrust


# Skip all tests if gzip support is not enabled
pytestmark = pytest.mark.skipif(
    not sdfrust.gzip_enabled(),
    reason="gzip feature not enabled"
)


# Test data
SIMPLE_SDF = """methane
  test    3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"""

SIMPLE_MOL2 = """@<TRIPOS>MOLECULE
methane
 5 4 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 MOL       0.0000
      2 H1          0.6289    0.6289    0.6289 H         1 MOL       0.0000
      3 H2         -0.6289   -0.6289    0.6289 H         1 MOL       0.0000
      4 H3         -0.6289    0.6289   -0.6289 H         1 MOL       0.0000
      5 H4          0.6289   -0.6289   -0.6289 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
"""

SIMPLE_XYZ = """3
water molecule
O  0.000000  0.000000  0.117300
H  0.756950  0.000000 -0.469200
H -0.756950  0.000000 -0.469200
"""


def create_gzip_file(content: str, suffix: str) -> str:
    """Create a gzip-compressed temporary file and return its path."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    with gzip.open(path, 'wt') as f:
        f.write(content)
    return path


class TestGzipEnabled:
    """Test gzip_enabled function."""

    def test_gzip_enabled_returns_true(self):
        """Test that gzip_enabled returns True when feature is enabled."""
        # This test only runs if gzip is enabled, so it should be True
        assert sdfrust.gzip_enabled() is True


class TestGzipSdfParsing:
    """Test parsing gzipped SDF files."""

    def test_parse_gzipped_sdf_file(self):
        """Test parsing a single molecule from a gzipped SDF file."""
        path = create_gzip_file(SIMPLE_SDF, ".sdf.gz")
        try:
            mol = sdfrust.parse_sdf_file(path)
            assert mol.name == "methane"
            assert mol.num_atoms == 5
            assert mol.num_bonds == 4
            assert mol.formula() == "CH4"
        finally:
            os.unlink(path)

    def test_parse_gzipped_sdf_file_multi(self):
        """Test parsing multiple molecules from a gzipped SDF file."""
        multi_sdf = SIMPLE_SDF + SIMPLE_SDF
        path = create_gzip_file(multi_sdf, ".sdf.gz")
        try:
            mols = sdfrust.parse_sdf_file_multi(path)
            assert len(mols) == 2
            assert mols[0].name == "methane"
            assert mols[1].name == "methane"
        finally:
            os.unlink(path)

    def test_iter_gzipped_sdf_file(self):
        """Test iterating over a gzipped SDF file."""
        multi_sdf = SIMPLE_SDF + SIMPLE_SDF + SIMPLE_SDF
        path = create_gzip_file(multi_sdf, ".sdf.gz")
        try:
            mols = list(sdfrust.iter_sdf_file(path))
            assert len(mols) == 3
            for mol in mols:
                assert mol.name == "methane"
        finally:
            os.unlink(path)


class TestGzipMol2Parsing:
    """Test parsing gzipped MOL2 files."""

    def test_parse_gzipped_mol2_file(self):
        """Test parsing a single molecule from a gzipped MOL2 file."""
        path = create_gzip_file(SIMPLE_MOL2, ".mol2.gz")
        try:
            mol = sdfrust.parse_mol2_file(path)
            assert mol.name == "methane"
            assert mol.num_atoms == 5
            assert mol.num_bonds == 4
        finally:
            os.unlink(path)

    def test_parse_gzipped_mol2_file_multi(self):
        """Test parsing multiple molecules from a gzipped MOL2 file."""
        multi_mol2 = SIMPLE_MOL2 + SIMPLE_MOL2
        path = create_gzip_file(multi_mol2, ".mol2.gz")
        try:
            mols = sdfrust.parse_mol2_file_multi(path)
            assert len(mols) == 2
        finally:
            os.unlink(path)

    def test_iter_gzipped_mol2_file(self):
        """Test iterating over a gzipped MOL2 file."""
        multi_mol2 = SIMPLE_MOL2 + SIMPLE_MOL2 + SIMPLE_MOL2
        path = create_gzip_file(multi_mol2, ".mol2.gz")
        try:
            mols = list(sdfrust.iter_mol2_file(path))
            assert len(mols) == 3
        finally:
            os.unlink(path)


class TestGzipXyzParsing:
    """Test parsing gzipped XYZ files."""

    def test_parse_gzipped_xyz_file(self):
        """Test parsing a single molecule from a gzipped XYZ file."""
        path = create_gzip_file(SIMPLE_XYZ, ".xyz.gz")
        try:
            mol = sdfrust.parse_xyz_file(path)
            assert mol.name == "water molecule"
            assert mol.num_atoms == 3
            assert mol.num_bonds == 0  # XYZ has no bonds
            assert mol.formula() == "H2O"
        finally:
            os.unlink(path)

    def test_parse_gzipped_xyz_file_multi(self):
        """Test parsing multiple molecules from a gzipped XYZ file."""
        multi_xyz = SIMPLE_XYZ + SIMPLE_XYZ
        path = create_gzip_file(multi_xyz, ".xyz.gz")
        try:
            mols = sdfrust.parse_xyz_file_multi(path)
            assert len(mols) == 2
        finally:
            os.unlink(path)

    def test_iter_gzipped_xyz_file(self):
        """Test iterating over a gzipped XYZ file."""
        multi_xyz = SIMPLE_XYZ + SIMPLE_XYZ + SIMPLE_XYZ
        path = create_gzip_file(multi_xyz, ".xyz.gz")
        try:
            mols = list(sdfrust.iter_xyz_file(path))
            assert len(mols) == 3
        finally:
            os.unlink(path)


class TestGzipAutoDetection:
    """Test automatic format detection with gzipped files."""

    def test_parse_auto_gzipped_sdf(self):
        """Test auto-detection parsing of gzipped SDF file."""
        path = create_gzip_file(SIMPLE_SDF, ".sdf.gz")
        try:
            mol = sdfrust.parse_auto_file(path)
            assert mol.name == "methane"
            assert mol.formula() == "CH4"
        finally:
            os.unlink(path)

    def test_parse_auto_gzipped_mol2(self):
        """Test auto-detection parsing of gzipped MOL2 file."""
        path = create_gzip_file(SIMPLE_MOL2, ".mol2.gz")
        try:
            mol = sdfrust.parse_auto_file(path)
            assert mol.name == "methane"
        finally:
            os.unlink(path)

    def test_parse_auto_gzipped_xyz(self):
        """Test auto-detection parsing of gzipped XYZ file."""
        path = create_gzip_file(SIMPLE_XYZ, ".xyz.gz")
        try:
            mol = sdfrust.parse_auto_file(path)
            assert mol.name == "water molecule"
            assert mol.formula() == "H2O"
        finally:
            os.unlink(path)

    def test_parse_auto_file_multi_gzipped(self):
        """Test auto-detection parsing of multiple molecules from gzipped file."""
        multi_sdf = SIMPLE_SDF + SIMPLE_SDF
        path = create_gzip_file(multi_sdf, ".sdf.gz")
        try:
            mols = sdfrust.parse_auto_file_multi(path)
            assert len(mols) == 2
        finally:
            os.unlink(path)


class TestGzipEdgeCases:
    """Test edge cases for gzip support."""

    def test_gzip_case_insensitive_extension(self):
        """Test that .GZ extension works (case insensitive)."""
        path = create_gzip_file(SIMPLE_SDF, ".sdf.GZ")
        try:
            mol = sdfrust.parse_sdf_file(path)
            assert mol.name == "methane"
        finally:
            os.unlink(path)

    def test_plain_file_still_works(self):
        """Test that plain (non-gzipped) files still work."""
        fd, path = tempfile.mkstemp(suffix=".sdf")
        try:
            os.write(fd, SIMPLE_SDF.encode())
            os.close(fd)
            mol = sdfrust.parse_sdf_file(path)
            assert mol.name == "methane"
        finally:
            os.unlink(path)
