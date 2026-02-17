"""Basic tests for sdfrust Python bindings."""

import os
import tempfile

import pytest

import sdfrust


# Path to test data (relative to workspace root)
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")


class TestVersion:
    """Test version information."""

    def test_version(self):
        """Test that version is accessible."""
        assert hasattr(sdfrust, "__version__")
        assert isinstance(sdfrust.__version__, str)
        assert sdfrust.__version__ == "0.5.0"


class TestAtom:
    """Test Atom class."""

    def test_create_atom(self):
        """Test creating an atom."""
        atom = sdfrust.Atom(0, "C", 1.0, 2.0, 3.0)
        assert atom.index == 0
        assert atom.element == "C"
        assert atom.x == 1.0
        assert atom.y == 2.0
        assert atom.z == 3.0
        assert atom.formal_charge == 0

    def test_atom_coords(self):
        """Test atom coordinate tuple."""
        atom = sdfrust.Atom(0, "N", 1.5, 2.5, 3.5)
        assert atom.coords() == (1.5, 2.5, 3.5)

    def test_atom_distance(self):
        """Test distance calculation between atoms."""
        atom1 = sdfrust.Atom(0, "C", 0.0, 0.0, 0.0)
        atom2 = sdfrust.Atom(1, "C", 3.0, 4.0, 0.0)
        assert abs(atom1.distance_to(atom2) - 5.0) < 1e-10

    def test_atom_charged(self):
        """Test is_charged method."""
        atom = sdfrust.Atom(0, "N", 0.0, 0.0, 0.0)
        assert not atom.is_charged()
        atom.formal_charge = 1
        assert atom.is_charged()

    def test_atom_repr(self):
        """Test atom string representation."""
        atom = sdfrust.Atom(0, "C", 1.0, 2.0, 3.0)
        assert "Atom" in repr(atom)
        assert "C" in repr(atom)


class TestBond:
    """Test Bond class."""

    def test_create_bond(self):
        """Test creating a bond."""
        order = sdfrust.BondOrder.single()
        bond = sdfrust.Bond(0, 1, order)
        assert bond.atom1 == 0
        assert bond.atom2 == 1

    def test_bond_order(self):
        """Test bond order types."""
        assert sdfrust.BondOrder.single().order() == 1.0
        assert sdfrust.BondOrder.double().order() == 2.0
        assert sdfrust.BondOrder.triple().order() == 3.0
        assert sdfrust.BondOrder.aromatic().order() == 1.5

    def test_bond_stereo(self):
        """Test bond stereochemistry."""
        stereo = sdfrust.BondStereo.up()
        assert str(stereo) == "up"

    def test_bond_contains_atom(self):
        """Test bond contains_atom method."""
        bond = sdfrust.Bond(0, 1, sdfrust.BondOrder.single())
        assert bond.contains_atom(0)
        assert bond.contains_atom(1)
        assert not bond.contains_atom(2)

    def test_bond_other_atom(self):
        """Test bond other_atom method."""
        bond = sdfrust.Bond(0, 1, sdfrust.BondOrder.single())
        assert bond.other_atom(0) == 1
        assert bond.other_atom(1) == 0
        assert bond.other_atom(2) is None


class TestMolecule:
    """Test Molecule class."""

    def test_create_molecule(self):
        """Test creating an empty molecule."""
        mol = sdfrust.Molecule("test")
        assert mol.name == "test"
        assert mol.num_atoms == 0
        assert mol.num_bonds == 0
        assert mol.is_empty()

    def test_add_atoms_bonds(self):
        """Test adding atoms and bonds."""
        mol = sdfrust.Molecule("water")
        mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))
        mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
        mol.add_bond(sdfrust.Bond(0, 2, sdfrust.BondOrder.single()))

        assert mol.num_atoms == 3
        assert mol.num_bonds == 2
        assert not mol.is_empty()

    def test_formula(self):
        """Test molecular formula calculation."""
        mol = sdfrust.Molecule("water")
        mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))
        assert mol.formula() == "H2O"

    def test_properties(self):
        """Test molecule properties."""
        mol = sdfrust.Molecule("test")
        mol.set_property("KEY", "value")
        assert mol.get_property("KEY") == "value"
        assert mol.get_property("NONEXISTENT") is None

    def test_neighbors(self):
        """Test finding atom neighbors."""
        mol = sdfrust.Molecule("methane")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "H", 1.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "H", -1.0, 0.0, 0.0))
        mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
        mol.add_bond(sdfrust.Bond(0, 2, sdfrust.BondOrder.single()))

        neighbors = mol.neighbors(0)
        assert 1 in neighbors
        assert 2 in neighbors
        assert len(neighbors) == 2


class TestSdfParsing:
    """Test SDF parsing functionality."""

    def test_parse_sdf_string(self):
        """Test parsing SDF from string."""
        sdf_content = (
            "methane\n"
            "\n"
            "\n"
            "  5  4  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    0.6289    0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -0.6289   -0.6289    0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "   -0.6289    0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    0.6289   -0.6289   -0.6289 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0  0  0  0\n"
            "  1  3  1  0  0  0  0\n"
            "  1  4  1  0  0  0  0\n"
            "  1  5  1  0  0  0  0\n"
            "M  END\n"
            "$$$$\n"
        )
        mol = sdfrust.parse_sdf_string(sdf_content)
        assert mol.name == "methane"
        assert mol.num_atoms == 5
        assert mol.num_bonds == 4
        assert mol.formula() == "CH4"

    def test_parse_sdf_file(self):
        """Test parsing SDF from file."""
        aspirin_path = os.path.join(TEST_DATA_DIR, "aspirin.sdf")
        if os.path.exists(aspirin_path):
            mol = sdfrust.parse_sdf_file(aspirin_path)
            assert mol.num_atoms > 0
            assert mol.num_bonds > 0

    def test_parse_sdf_file_multi(self):
        """Test parsing multiple molecules from SDF file."""
        caffeine_path = os.path.join(TEST_DATA_DIR, "caffeine.sdf")
        if os.path.exists(caffeine_path):
            mols = sdfrust.parse_sdf_file_multi(caffeine_path)
            assert len(mols) >= 1

    def test_invalid_sdf(self):
        """Test that invalid SDF raises an error."""
        with pytest.raises(ValueError):
            sdfrust.parse_sdf_string("invalid content")


class TestMol2Parsing:
    """Test MOL2 parsing functionality."""

    def test_parse_mol2_string(self):
        """Test parsing MOL2 from string."""
        mol2_content = """@<TRIPOS>MOLECULE
water
 3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
      1 O1          0.0000    0.0000    0.0000 O.3       1 MOL       0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 MOL       0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 MOL       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
"""
        mol = sdfrust.parse_mol2_string(mol2_content)
        assert mol.name == "water"
        assert mol.num_atoms == 3
        assert mol.num_bonds == 2
        assert mol.formula() == "H2O"

    def test_parse_mol2_file(self):
        """Test parsing MOL2 from file."""
        benzene_path = os.path.join(TEST_DATA_DIR, "benzene.mol2")
        if os.path.exists(benzene_path):
            mol = sdfrust.parse_mol2_file(benzene_path)
            assert mol.num_atoms > 0


class TestSdfWriting:
    """Test SDF writing functionality."""

    def test_write_sdf_string(self):
        """Test writing SDF to string."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "O", 1.2, 0.0, 0.0))
        mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.double()))

        sdf_string = sdfrust.write_sdf_string(mol)
        assert "test" in sdf_string
        assert "V2000" in sdf_string
        assert "$$$$" in sdf_string

    def test_write_sdf_file(self):
        """Test writing SDF to file."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
            temp_path = f.name

        try:
            sdfrust.write_sdf_file(mol, temp_path)
            # Read back and verify
            mol2 = sdfrust.parse_sdf_file(temp_path)
            assert mol2.name == "test"
            assert mol2.num_atoms == 1
        finally:
            os.unlink(temp_path)


class TestXyzParsing:
    """Test XYZ parsing functionality."""

    def test_parse_xyz_string(self):
        """Test parsing XYZ from string."""
        xyz_content = """3
water molecule
O  0.000000  0.000000  0.117300
H  0.756950  0.000000 -0.469200
H -0.756950  0.000000 -0.469200
"""
        mol = sdfrust.parse_xyz_string(xyz_content)
        assert mol.name == "water molecule"
        assert mol.num_atoms == 3
        assert mol.num_bonds == 0  # XYZ has no bonds
        assert mol.formula() == "H2O"

    def test_parse_xyz_file(self):
        """Test parsing XYZ from file."""
        water_path = os.path.join(TEST_DATA_DIR, "water.xyz")
        if os.path.exists(water_path):
            mol = sdfrust.parse_xyz_file(water_path)
            assert mol.name == "water molecule"
            assert mol.num_atoms == 3
            assert mol.formula() == "H2O"

    def test_parse_xyz_file_multi(self):
        """Test parsing multiple molecules from XYZ file."""
        multi_path = os.path.join(TEST_DATA_DIR, "multi.xyz")
        if os.path.exists(multi_path):
            mols = sdfrust.parse_xyz_file_multi(multi_path)
            assert len(mols) == 3
            assert mols[0].name == "water molecule"
            assert mols[1].name == "methane"
            assert mols[2].name == "diatomic hydrogen"

    def test_parse_xyz_string_multi(self):
        """Test parsing multiple molecules from XYZ string."""
        xyz_content = """2
mol1
C  0.0  0.0  0.0
O  1.2  0.0  0.0
3
mol2
N  0.0  0.0  0.0
H  1.0  0.0  0.0
H -1.0  0.0  0.0
"""
        mols = sdfrust.parse_xyz_string_multi(xyz_content)
        assert len(mols) == 2
        assert mols[0].name == "mol1"
        assert mols[1].name == "mol2"

    def test_iter_xyz_file(self):
        """Test iterating over XYZ file."""
        multi_path = os.path.join(TEST_DATA_DIR, "multi.xyz")
        if os.path.exists(multi_path):
            count = 0
            for mol in sdfrust.iter_xyz_file(multi_path):
                count += 1
                assert mol.num_atoms > 0
            assert count == 3

    def test_invalid_xyz(self):
        """Test that invalid XYZ raises an error."""
        with pytest.raises(ValueError):
            sdfrust.parse_xyz_string("invalid content")

    def test_detect_format_xyz(self):
        """Test format detection for XYZ content."""
        xyz_content = "3\nwater\nO 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n"
        fmt = sdfrust.detect_format(xyz_content)
        assert fmt == "xyz"

    def test_parse_auto_xyz(self):
        """Test auto-detection parsing of XYZ content."""
        xyz_content = """3
water
O  0.0  0.0  0.0
H  1.0  0.0  0.0
H -1.0  0.0  0.0
"""
        mol = sdfrust.parse_auto_string(xyz_content)
        assert mol.name == "water"
        assert mol.num_atoms == 3


class TestIterators:
    """Test iterator functionality."""

    def test_iter_sdf_file(self):
        """Test iterating over SDF file."""
        aspirin_path = os.path.join(TEST_DATA_DIR, "aspirin.sdf")
        if os.path.exists(aspirin_path):
            count = 0
            for mol in sdfrust.iter_sdf_file(aspirin_path):
                count += 1
                assert mol.num_atoms > 0
            assert count >= 1


class TestDescriptors:
    """Test molecular descriptor calculations."""

    def test_molecular_weight(self):
        """Test molecular weight calculation."""
        mol = sdfrust.Molecule("water")
        mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))

        mw = mol.molecular_weight()
        assert mw is not None
        assert abs(mw - 18.015) < 0.01

    def test_heavy_atom_count(self):
        """Test heavy atom count."""
        mol = sdfrust.Molecule("water")
        mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))

        assert mol.heavy_atom_count() == 1  # Only oxygen


class TestGeometry:
    """Test geometry operations."""

    def test_centroid(self):
        """Test centroid calculation."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "C", 2.0, 0.0, 0.0))

        centroid = mol.centroid()
        assert centroid is not None
        x, y, z = centroid
        assert abs(x - 1.0) < 1e-10
        assert abs(y - 0.0) < 1e-10
        assert abs(z - 0.0) < 1e-10

    def test_translate(self):
        """Test molecule translation."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))

        mol.translate(1.0, 2.0, 3.0)
        atom = mol.get_atom(0)
        assert abs(atom.x - 1.0) < 1e-10
        assert abs(atom.y - 2.0) < 1e-10
        assert abs(atom.z - 3.0) < 1e-10

    def test_center(self):
        """Test centering molecule at origin."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 2.0, 2.0, 2.0))
        mol.add_atom(sdfrust.Atom(1, "C", 4.0, 4.0, 4.0))

        mol.center()
        centroid = mol.centroid()
        assert centroid is not None
        x, y, z = centroid
        assert abs(x) < 1e-10
        assert abs(y) < 1e-10
        assert abs(z) < 1e-10


# NumPy tests (only run if numpy is available)
try:
    import numpy as np

    class TestNumpy:
        """Test NumPy integration."""

        def test_get_coords_array(self):
            """Test getting coordinates as NumPy array."""
            mol = sdfrust.Molecule("test")
            mol.add_atom(sdfrust.Atom(0, "C", 1.0, 2.0, 3.0))
            mol.add_atom(sdfrust.Atom(1, "O", 4.0, 5.0, 6.0))

            coords = mol.get_coords_array()
            assert coords.shape == (2, 3)
            assert np.allclose(coords[0], [1.0, 2.0, 3.0])
            assert np.allclose(coords[1], [4.0, 5.0, 6.0])

        def test_set_coords_array(self):
            """Test setting coordinates from NumPy array."""
            mol = sdfrust.Molecule("test")
            mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
            mol.add_atom(sdfrust.Atom(1, "O", 0.0, 0.0, 0.0))

            new_coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
            mol.set_coords_array(new_coords)

            atom0 = mol.get_atom(0)
            atom1 = mol.get_atom(1)
            assert abs(atom0.x - 1.0) < 1e-10
            assert abs(atom1.z - 6.0) < 1e-10

        def test_get_atomic_numbers(self):
            """Test getting atomic numbers."""
            mol = sdfrust.Molecule("water")
            mol.add_atom(sdfrust.Atom(0, "O", 0.0, 0.0, 0.0))
            mol.add_atom(sdfrust.Atom(1, "H", 0.96, 0.0, 0.0))
            mol.add_atom(sdfrust.Atom(2, "H", -0.24, 0.93, 0.0))

            atomic_nums = mol.get_atomic_numbers()
            assert list(atomic_nums) == [8, 1, 1]  # O, H, H

except ImportError:
    pass  # NumPy not available, skip these tests


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
