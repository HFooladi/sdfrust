"""Tests for geometry operations in sdfrust Python bindings.

These tests require the geometry feature to be enabled:
    maturin develop --features numpy,geometry
"""

import math

import pytest

import sdfrust


def has_geometry_feature():
    """Check if geometry feature is enabled."""
    mol = sdfrust.Molecule("test")
    return hasattr(mol, "distance_matrix")


# Skip all tests if geometry feature is not enabled
pytestmark = pytest.mark.skipif(
    not has_geometry_feature(),
    reason="geometry feature not enabled",
)


class TestDistanceMatrix:
    """Test distance matrix calculation."""

    def test_distance_matrix_empty(self):
        """Test distance matrix of empty molecule."""
        mol = sdfrust.Molecule("empty")
        matrix = mol.distance_matrix()
        assert matrix == []

    def test_distance_matrix_single_atom(self):
        """Test distance matrix with single atom."""
        mol = sdfrust.Molecule("single")
        mol.add_atom(sdfrust.Atom(0, "C", 1.0, 2.0, 3.0))

        matrix = mol.distance_matrix()
        assert len(matrix) == 1
        assert len(matrix[0]) == 1
        assert abs(matrix[0][0]) < 1e-10

    def test_distance_matrix_two_atoms(self):
        """Test distance matrix with two atoms."""
        mol = sdfrust.Molecule("two")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "C", 3.0, 4.0, 0.0))  # Distance = 5

        matrix = mol.distance_matrix()
        assert len(matrix) == 2
        assert abs(matrix[0][1] - 5.0) < 1e-10
        assert abs(matrix[1][0] - 5.0) < 1e-10
        assert abs(matrix[0][0]) < 1e-10
        assert abs(matrix[1][1]) < 1e-10

    def test_distance_matrix_symmetric(self):
        """Test that distance matrix is symmetric."""
        mol = sdfrust.Molecule("triangle")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(2, "C", 0.5, 0.866, 0.0))

        matrix = mol.distance_matrix()

        for i in range(3):
            for j in range(3):
                assert abs(matrix[i][j] - matrix[j][i]) < 1e-10


class TestRMSD:
    """Test RMSD calculation."""

    def test_rmsd_identical_molecules(self):
        """Test RMSD of identical molecules is zero."""
        mol1 = sdfrust.Molecule("mol1")
        mol1.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol1.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))

        mol2 = sdfrust.Molecule("mol2")
        mol2.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol2.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))

        rmsd = mol1.rmsd_to(mol2)
        assert abs(rmsd) < 1e-10

    def test_rmsd_translated_molecule(self):
        """Test RMSD of translated molecule."""
        mol1 = sdfrust.Molecule("mol1")
        mol1.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol1.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))

        mol2 = sdfrust.Molecule("mol2")
        mol2.add_atom(sdfrust.Atom(0, "C", 1.0, 0.0, 0.0))  # Translated by 1
        mol2.add_atom(sdfrust.Atom(1, "C", 2.0, 0.0, 0.0))

        rmsd = mol1.rmsd_to(mol2)
        assert abs(rmsd - 1.0) < 1e-10

    def test_rmsd_different_atom_counts(self):
        """Test RMSD fails with different atom counts."""
        mol1 = sdfrust.Molecule("mol1")
        mol1.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))

        mol2 = sdfrust.Molecule("mol2")
        mol2.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol2.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))

        with pytest.raises(ValueError):
            mol1.rmsd_to(mol2)


class TestRotation:
    """Test rotation operations."""

    def test_rotate_identity(self):
        """Test rotation by zero angle."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 1.0, 2.0, 3.0))

        mol.rotate([0.0, 0.0, 1.0], 0.0)

        atom = mol.get_atom(0)
        assert abs(atom.x - 1.0) < 1e-10
        assert abs(atom.y - 2.0) < 1e-10
        assert abs(atom.z - 3.0) < 1e-10

    def test_rotate_90_degrees_z_axis(self):
        """Test 90 degree rotation around Z axis."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 1.0, 0.0, 0.0))

        mol.rotate([0.0, 0.0, 1.0], math.pi / 2)

        # (1, 0, 0) -> (0, 1, 0)
        atom = mol.get_atom(0)
        assert abs(atom.x) < 1e-10
        assert abs(atom.y - 1.0) < 1e-10
        assert abs(atom.z) < 1e-10

    def test_rotate_preserves_distances(self):
        """Test that rotation preserves interatomic distances."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
        mol.add_atom(sdfrust.Atom(1, "C", 1.0, 0.0, 0.0))

        distance_before = mol.get_atom(0).distance_to(mol.get_atom(1))

        mol.rotate([1.0, 2.0, 3.0], 1.234)

        distance_after = mol.get_atom(0).distance_to(mol.get_atom(1))

        assert abs(distance_before - distance_after) < 1e-10


class TestTransform:
    """Test transform operations."""

    def test_apply_rotation_matrix_identity(self):
        """Test applying identity rotation matrix."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 1.0, 2.0, 3.0))

        identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        mol.apply_rotation_matrix(identity)

        atom = mol.get_atom(0)
        assert abs(atom.x - 1.0) < 1e-10
        assert abs(atom.y - 2.0) < 1e-10
        assert abs(atom.z - 3.0) < 1e-10

    def test_apply_transform(self):
        """Test applying rotation + translation."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 1.0, 0.0, 0.0))

        # 90° rotation around Z, then translate by (0, 0, 5)
        rot_90_z = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        mol.apply_transform(rot_90_z, [0.0, 0.0, 5.0])

        # (1, 0, 0) -> (0, 1, 0) -> (0, 1, 5)
        atom = mol.get_atom(0)
        assert abs(atom.x) < 1e-10
        assert abs(atom.y - 1.0) < 1e-10
        assert abs(atom.z - 5.0) < 1e-10


class TestIntegration:
    """Integration tests combining multiple operations."""

    def test_center_and_rotate(self):
        """Test centering then rotating a molecule."""
        mol = sdfrust.Molecule("test")
        mol.add_atom(sdfrust.Atom(0, "C", 10.0, 10.0, 10.0))
        mol.add_atom(sdfrust.Atom(1, "C", 11.0, 10.0, 10.0))

        # Center first
        mol.center()
        centroid = mol.centroid()
        assert centroid is not None
        assert abs(centroid[0]) < 1e-10
        assert abs(centroid[1]) < 1e-10
        assert abs(centroid[2]) < 1e-10

        # Rotate
        mol.rotate([0.0, 0.0, 1.0], math.pi / 2)

        # Centroid should still be at origin
        centroid = mol.centroid()
        assert abs(centroid[0]) < 1e-10
        assert abs(centroid[1]) < 1e-10
        assert abs(centroid[2]) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
