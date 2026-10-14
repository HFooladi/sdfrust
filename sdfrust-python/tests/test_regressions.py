"""Python-level regression tests for bugs found in the v0.6.0 audit."""

import os

import pytest

import sdfrust

TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")


def _read(name):
    with open(os.path.join(TEST_DATA, name)) as f:
        return f.read()


def test_v3000_through_parse_sdf_string():
    mol = sdfrust.parse_sdf_string(_read("v3000_benzene.sdf"))
    assert mol.num_atoms == 6


def test_trailing_blank_lines():
    assert len(sdfrust.parse_sdf_string_multi(_read("aspirin.sdf") + "\n\n")) == 1


def test_isotope_and_radical_attributes():
    block = (
        "x\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\n"
        "M  ISO  1   1 235\nM  RAD  1   1   2\nM  END\n$$$$\n"
    )
    atom = sdfrust.parse_sdf_string(block).atoms[0]
    assert atom.isotope == 235
    assert atom.radical == 2
    assert atom.radical_electrons == 1


def test_add_bond_out_of_range_raises_index_error():
    mol = sdfrust.Molecule("x")
    mol.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
    with pytest.raises(IndexError):
        mol.add_bond(sdfrust.Bond(0, 7, sdfrust.BondOrder.single()))


def test_ecfp_zero_bits_raises_value_error():
    mol = sdfrust.parse_sdf_string(_read("aspirin.sdf"))
    with pytest.raises(ValueError):
        mol.ecfp(2, 0)


def test_v2000_writer_refuses_too_many_atoms():
    mol = sdfrust.Molecule("big")
    for i in range(1000):
        mol.add_atom(sdfrust.Atom(i, "C", float(i % 100), 0.0, 0.0))
    with pytest.raises(ValueError):
        sdfrust.write_sdf_string(mol)
    assert "V3000" in sdfrust.write_sdf_auto_string(mol)


def test_ogb_atomic_number_is_z_minus_one():
    mol = sdfrust.parse_sdf_string(_read("aspirin.sdf"))
    feats = mol.ogb_atom_features()
    for atom, f in zip(mol.atoms, feats):
        expected = {"C": 5, "O": 7, "H": 0}[atom.element]
        assert f[0] == expected
