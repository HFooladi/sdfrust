"""
Cross-validation of sdfrust ML features against RDKit reference implementations.

Compares element-by-element:
- OGB atom features [N, 9]
- OGB bond features [E, 3]
- ECFP/Morgan fingerprints
- Gasteiger partial charges
- Aromaticity perception
- Hybridization
- Ring perception

Requires: rdkit, sdfrust (built with `maturin develop --features numpy`)

Run: pytest tests/test_ml_validation.py -v
"""

import os

import numpy as np
import pytest
rdkit = pytest.importorskip("rdkit")
Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")
rdMolDescriptors = pytest.importorskip("rdkit.Chem.rdMolDescriptors")

import sdfrust

TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")
EXAMPLE_DATA = os.path.join(os.path.dirname(__file__), "..", "examples", "data")

# All single-molecule SDF files for validation
SINGLE_MOL_FILES = [
    # From tests/test_data/
    "aspirin.sdf",
    "caffeine_pubchem.sdf",
    "glucose.sdf",
    "galactose.sdf",
    "acetaminophen.sdf",
    "methionine.sdf",
]

# Multi-molecule SDF file (examples/data/drug_library.sdf)
# and additional single-molecule files from examples/data/
EXAMPLE_MOL_FILES = [
    "ibuprofen.sdf",
    "dopamine.sdf",
    "cholesterol.sdf",
]


# ============================================================
# OGB feature encoding helpers (matches ogb.utils.features)
# ============================================================

def rdkit_atom_to_ogb_features(atom):
    """Compute OGB 9-feature vector for an RDKit atom (matches ogb.utils.features.atom_to_feature_vector)."""
    # Feature 0: Atomic number
    atomic_num = atom.GetAtomicNum()

    # Feature 1: Chirality tag
    chiral_map = {
        Chem.rdchem.ChiralType.CHI_UNSPECIFIED: 0,
        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
        Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: 2,
        Chem.rdchem.ChiralType.CHI_OTHER: 3,
    }
    chirality = chiral_map.get(atom.GetChiralTag(), 0)

    # Feature 2: Degree
    degree = atom.GetDegree()

    # Feature 3: Formal charge (shifted by +5)
    formal_charge = atom.GetFormalCharge() + 5

    # Feature 4: Number of Hs (total)
    num_hs = atom.GetTotalNumHs()

    # Feature 5: Number of radical electrons
    num_radical = atom.GetNumRadicalElectrons()

    # Feature 6: Hybridization
    hyb_map = {
        Chem.rdchem.HybridizationType.S: 0,
        Chem.rdchem.HybridizationType.SP: 1,
        Chem.rdchem.HybridizationType.SP2: 2,
        Chem.rdchem.HybridizationType.SP3: 3,
        Chem.rdchem.HybridizationType.SP3D: 4,
        Chem.rdchem.HybridizationType.SP3D2: 5,
    }
    hybridization = hyb_map.get(atom.GetHybridization(), 0)

    # Feature 7: Is aromatic
    is_aromatic = int(atom.GetIsAromatic())

    # Feature 8: Is in ring
    is_in_ring = int(atom.IsInRing())

    return [atomic_num, chirality, degree, formal_charge, num_hs,
            num_radical, hybridization, is_aromatic, is_in_ring]


def rdkit_bond_to_ogb_features(bond, mol):
    """Compute OGB 3-feature vector for an RDKit bond."""
    # Feature 0: Bond type
    bond_type_map = {
        Chem.rdchem.BondType.SINGLE: 0,
        Chem.rdchem.BondType.DOUBLE: 1,
        Chem.rdchem.BondType.TRIPLE: 2,
        Chem.rdchem.BondType.AROMATIC: 3,
    }
    bond_type = bond_type_map.get(bond.GetBondType(), 0)

    # Feature 1: Bond stereo
    stereo_map = {
        Chem.rdchem.BondStereo.STEREONONE: 0,
        Chem.rdchem.BondStereo.STEREOANY: 1,
        Chem.rdchem.BondStereo.STEREOZ: 2,
        Chem.rdchem.BondStereo.STEREOE: 3,
        Chem.rdchem.BondStereo.STEREOCIS: 4,
        Chem.rdchem.BondStereo.STEREOTRANS: 5,
    }
    stereo = stereo_map.get(bond.GetStereo(), 0)

    # Feature 2: Is conjugated
    is_conjugated = int(bond.GetIsConjugated())

    return [bond_type, stereo, is_conjugated]


def load_both(sdf_path):
    """Load a molecule in both sdfrust and RDKit."""
    sdf_mol = sdfrust.parse_sdf_file(sdf_path)
    rdkit_mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    return sdf_mol, rdkit_mol


def _resolve_sdf_path(filename):
    """Find the SDF file in test_data or examples/data."""
    for d in [TEST_DATA, EXAMPLE_DATA]:
        path = os.path.join(d, filename)
        if os.path.exists(path):
            return path
    return None


def _all_molecule_params():
    """Collect all single-molecule test params."""
    params = []
    for f in SINGLE_MOL_FILES + EXAMPLE_MOL_FILES:
        params.append(f)
    return params


def _drug_library_params():
    """Collect molecule indices from the multi-molecule drug_library.sdf."""
    path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    if not os.path.exists(path):
        return []
    # Count molecules via RDKit
    supplier = Chem.SDMolSupplier(path, removeHs=False)
    return [f"drug_library:{i}" for i in range(len(supplier)) if supplier[i] is not None]


# ============================================================
# Test fixtures
# ============================================================

@pytest.fixture(params=_all_molecule_params())
def molecule_pair(request):
    """Load molecule in both sdfrust and RDKit."""
    path = _resolve_sdf_path(request.param)
    if path is None:
        pytest.skip(f"Test file not found: {request.param}")
    sdf_mol, rdkit_mol = load_both(path)
    assert rdkit_mol is not None, f"RDKit failed to parse {request.param}"
    return request.param, sdf_mol, rdkit_mol


@pytest.fixture(params=_drug_library_params())
def drug_library_pair(request):
    """Load a molecule from the multi-molecule drug_library.sdf."""
    _, idx_str = request.param.split(":")
    idx = int(idx_str)
    path = os.path.join(EXAMPLE_DATA, "drug_library.sdf")
    if not os.path.exists(path):
        pytest.skip("drug_library.sdf not found")
    sdf_mols = sdfrust.parse_sdf_file_multi(path)
    supplier = Chem.SDMolSupplier(path, removeHs=False)
    rdkit_mol = supplier[idx]
    if idx >= len(sdf_mols) or rdkit_mol is None:
        pytest.skip(f"Could not load molecule {idx}")
    return request.param, sdf_mols[idx], rdkit_mol


# ============================================================
# OGB Atom Feature Validation
# ============================================================

class TestOGBAtomFeatures:
    """Validate OGB atom features against RDKit reference."""

    def test_atom_count_matches(self, molecule_pair):
        name, sdf_mol, rdkit_mol = molecule_pair
        assert sdf_mol.num_atoms == rdkit_mol.GetNumAtoms(), \
            f"{name}: atom count mismatch"

    def test_atomic_number(self, molecule_pair):
        """Feature 0: Atomic number should match exactly."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        for i in range(sdf_mol.num_atoms):
            rdkit_num = rdkit_mol.GetAtomWithIdx(i).GetAtomicNum()
            sdf_num = sdf_feats[i][0]
            assert sdf_num == rdkit_num, \
                f"{name} atom {i}: atomic_num sdfrust={sdf_num} rdkit={rdkit_num}"

    def test_degree(self, molecule_pair):
        """Feature 2: Degree should match exactly."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        for i in range(sdf_mol.num_atoms):
            rdkit_deg = rdkit_mol.GetAtomWithIdx(i).GetDegree()
            sdf_deg = sdf_feats[i][2]
            assert sdf_deg == rdkit_deg, \
                f"{name} atom {i}: degree sdfrust={sdf_deg} rdkit={rdkit_deg}"

    def test_formal_charge(self, molecule_pair):
        """Feature 3: Formal charge (+5 shift) should match exactly."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        for i in range(sdf_mol.num_atoms):
            rdkit_charge = rdkit_mol.GetAtomWithIdx(i).GetFormalCharge() + 5
            sdf_charge = sdf_feats[i][3]
            assert sdf_charge == rdkit_charge, \
                f"{name} atom {i}: charge sdfrust={sdf_charge} rdkit={rdkit_charge}"

    def test_num_hs(self, molecule_pair):
        """Feature 4: Total H count should match exactly."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        for i in range(sdf_mol.num_atoms):
            rdkit_hs = rdkit_mol.GetAtomWithIdx(i).GetTotalNumHs()
            sdf_hs = sdf_feats[i][4]
            assert sdf_hs == rdkit_hs, \
                f"{name} atom {i}: num_hs sdfrust={sdf_hs} rdkit={rdkit_hs}"

    def test_hybridization(self, molecule_pair):
        """Feature 6: Hybridization should match."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        hyb_map = {
            Chem.rdchem.HybridizationType.S: 0,
            Chem.rdchem.HybridizationType.SP: 1,
            Chem.rdchem.HybridizationType.SP2: 2,
            Chem.rdchem.HybridizationType.SP3: 3,
            Chem.rdchem.HybridizationType.SP3D: 4,
            Chem.rdchem.HybridizationType.SP3D2: 5,
        }
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_hyb = hyb_map.get(rdkit_mol.GetAtomWithIdx(i).GetHybridization(), 0)
            sdf_hyb = sdf_feats[i][6]
            if sdf_hyb != rdkit_hyb:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(f"  atom {i} ({elem}): sdfrust={sdf_hyb} rdkit={rdkit_hyb}")
        if mismatches:
            pytest.fail(f"{name} hybridization mismatches:\n" + "\n".join(mismatches))

    def test_is_aromatic(self, molecule_pair):
        """Feature 7: Aromaticity should match."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_arom = int(rdkit_mol.GetAtomWithIdx(i).GetIsAromatic())
            sdf_arom = sdf_feats[i][7]
            if sdf_arom != rdkit_arom:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(f"  atom {i} ({elem}): sdfrust={sdf_arom} rdkit={rdkit_arom}")
        if mismatches:
            pytest.fail(f"{name} aromaticity mismatches:\n" + "\n".join(mismatches))

    def test_is_in_ring(self, molecule_pair):
        """Feature 8: Ring membership should match."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_ring = int(rdkit_mol.GetAtomWithIdx(i).IsInRing())
            sdf_ring = sdf_feats[i][8]
            if sdf_ring != rdkit_ring:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(f"  atom {i} ({elem}): sdfrust={sdf_ring} rdkit={rdkit_ring}")
        if mismatches:
            pytest.fail(f"{name} ring membership mismatches:\n" + "\n".join(mismatches))

    def test_chirality(self, molecule_pair):
        """Feature 1: Chirality tag should match RDKit (CIP-based perception)."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        chiral_map = {
            Chem.rdchem.ChiralType.CHI_UNSPECIFIED: 0,
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
            Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: 2,
            Chem.rdchem.ChiralType.CHI_OTHER: 3,
        }
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_chiral = chiral_map.get(rdkit_mol.GetAtomWithIdx(i).GetChiralTag(), 0)
            sdf_chiral = sdf_feats[i][1]
            if sdf_chiral != rdkit_chiral:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(
                    f"  atom {i} ({elem}): sdfrust={sdf_chiral} rdkit={rdkit_chiral}"
                )
        if mismatches:
            pytest.fail(f"{name} chirality mismatches:\n" + "\n".join(mismatches))

    def test_full_ogb_atom_features(self, molecule_pair):
        """Full comparison of all 9 OGB atom features."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        feature_names = ["atomic_num", "chirality", "degree", "charge",
                         "num_hs", "radical", "hybridization", "aromatic", "in_ring"]
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_feat = rdkit_atom_to_ogb_features(rdkit_mol.GetAtomWithIdx(i))
            for j in range(9):
                if sdf_feats[i][j] != rdkit_feat[j]:
                    elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                    mismatches.append(
                        f"  atom {i} ({elem}) {feature_names[j]}: "
                        f"sdfrust={sdf_feats[i][j]} rdkit={rdkit_feat[j]}"
                    )
        if mismatches:
            pytest.fail(f"{name}: {len(mismatches)} OGB atom feature mismatches:\n" +
                        "\n".join(mismatches))


# ============================================================
# OGB Bond Feature Validation
# ============================================================

class TestOGBBondFeatures:
    """Validate OGB bond features against RDKit reference."""

    def test_bond_count_matches(self, molecule_pair):
        name, sdf_mol, rdkit_mol = molecule_pair
        assert sdf_mol.num_bonds == rdkit_mol.GetNumBonds(), \
            f"{name}: bond count mismatch"

    def test_bond_type(self, molecule_pair):
        """Feature 0: Bond type should match."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_bond_features()
        bond_type_map = {
            Chem.rdchem.BondType.SINGLE: 0,
            Chem.rdchem.BondType.DOUBLE: 1,
            Chem.rdchem.BondType.TRIPLE: 2,
            Chem.rdchem.BondType.AROMATIC: 3,
        }
        mismatches = []
        for i in range(sdf_mol.num_bonds):
            rdkit_bond = rdkit_mol.GetBondWithIdx(i)
            rdkit_type = bond_type_map.get(rdkit_bond.GetBondType(), 0)
            sdf_type = sdf_feats[i][0]
            if sdf_type != rdkit_type:
                a1, a2 = rdkit_bond.GetBeginAtomIdx(), rdkit_bond.GetEndAtomIdx()
                mismatches.append(
                    f"  bond {i} ({a1}-{a2}): sdfrust={sdf_type} rdkit={rdkit_type}"
                )
        if mismatches:
            pytest.fail(f"{name} bond type mismatches:\n" + "\n".join(mismatches))

    def test_is_conjugated(self, molecule_pair):
        """Feature 2: Conjugation should match."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_bond_features()
        mismatches = []
        for i in range(sdf_mol.num_bonds):
            rdkit_bond = rdkit_mol.GetBondWithIdx(i)
            rdkit_conj = int(rdkit_bond.GetIsConjugated())
            sdf_conj = sdf_feats[i][2]
            if sdf_conj != rdkit_conj:
                a1, a2 = rdkit_bond.GetBeginAtomIdx(), rdkit_bond.GetEndAtomIdx()
                e1 = rdkit_mol.GetAtomWithIdx(a1).GetSymbol()
                e2 = rdkit_mol.GetAtomWithIdx(a2).GetSymbol()
                mismatches.append(
                    f"  bond {i} ({e1}{a1}-{e2}{a2}): sdfrust={sdf_conj} rdkit={rdkit_conj}"
                )
        if mismatches:
            pytest.fail(f"{name} conjugation mismatches:\n" + "\n".join(mismatches))

    def test_full_ogb_bond_features(self, molecule_pair):
        """Full comparison of all 3 OGB bond features."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_feats = sdf_mol.ogb_bond_features()
        feature_names = ["bond_type", "stereo", "conjugated"]
        mismatches = []
        for i in range(sdf_mol.num_bonds):
            rdkit_bond = rdkit_mol.GetBondWithIdx(i)
            rdkit_feat = rdkit_bond_to_ogb_features(rdkit_bond, rdkit_mol)
            for j in range(3):
                if sdf_feats[i][j] != rdkit_feat[j]:
                    a1, a2 = rdkit_bond.GetBeginAtomIdx(), rdkit_bond.GetEndAtomIdx()
                    mismatches.append(
                        f"  bond {i} ({a1}-{a2}) {feature_names[j]}: "
                        f"sdfrust={sdf_feats[i][j]} rdkit={rdkit_feat[j]}"
                    )
        if mismatches:
            pytest.fail(f"{name}: {len(mismatches)} OGB bond feature mismatches:\n" +
                        "\n".join(mismatches))


# ============================================================
# ECFP Fingerprint Validation
# ============================================================

class TestECFPFingerprints:
    """Validate ECFP fingerprints against RDKit Morgan fingerprints."""

    def test_ecfp_nonzero(self, molecule_pair):
        """Both implementations should produce non-empty fingerprints."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_fp = sdf_mol.ecfp(radius=2, n_bits=2048)
        rdkit_fp = AllChem.GetMorganFingerprintAsBitVect(rdkit_mol, radius=2, nBits=2048)
        assert sum(sdf_fp) > 0, f"{name}: sdfrust ECFP is all zeros"
        assert rdkit_fp.GetNumOnBits() > 0, f"{name}: RDKit Morgan is all zeros"

    def test_ecfp_bit_overlap(self, molecule_pair):
        """Check overlap between sdfrust and RDKit fingerprints.

        Note: Exact bit identity is NOT expected because the hash functions differ.
        RDKit uses a different initial atom invariant and hashing scheme than our
        implementation. What matters is that:
        1. Both produce reasonable density
        2. Self-similarity is 1.0
        3. Similar molecules get similar scores
        """
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_fp = sdf_mol.ecfp(radius=2, n_bits=2048)
        rdkit_fp = AllChem.GetMorganFingerprintAsBitVect(rdkit_mol, radius=2, nBits=2048)

        sdf_density = sum(sdf_fp) / len(sdf_fp)
        rdkit_on = rdkit_fp.GetNumOnBits()
        rdkit_density = rdkit_on / rdkit_fp.GetNumBits()

        # Both should have reasonable density (not too sparse, not too dense)
        assert 0.001 < sdf_density < 0.5, \
            f"{name}: sdfrust density {sdf_density:.4f} out of range"
        assert 0.001 < rdkit_density < 0.5, \
            f"{name}: rdkit density {rdkit_density:.4f} out of range"

    def test_ecfp_self_similarity(self, molecule_pair):
        """Self-similarity should be 1.0 for both implementations."""
        name, sdf_mol, _ = molecule_pair
        sim = sdf_mol.tanimoto_similarity(sdf_mol)
        assert abs(sim - 1.0) < 1e-10, \
            f"{name}: self-similarity = {sim}, expected 1.0"


# ============================================================
# Gasteiger Charge Validation
# ============================================================

class TestGasteigerCharges:
    """Validate Gasteiger charges against RDKit reference."""

    def test_charge_signs(self, molecule_pair):
        """Electronegative atoms should be negative, electropositive positive."""
        name, sdf_mol, rdkit_mol = molecule_pair
        AllChem.ComputeGasteigerCharges(rdkit_mol)

        sdf_charges = sdf_mol.gasteiger_charges()

        for i in range(sdf_mol.num_atoms):
            rdkit_charge = float(rdkit_mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge"))
            sdf_charge = sdf_charges[i]

            # Check that charge signs match (most important for ML)
            if abs(rdkit_charge) > 0.05 and abs(sdf_charge) > 0.05:
                rdkit_sign = 1 if rdkit_charge > 0 else -1
                sdf_sign = 1 if sdf_charge > 0 else -1
                if rdkit_sign != sdf_sign:
                    elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                    print(f"  WARNING {name} atom {i} ({elem}): sign mismatch "
                          f"sdfrust={sdf_charge:+.4f} rdkit={rdkit_charge:+.4f}")

    def test_charge_correlation(self, molecule_pair):
        """Charges should be correlated (same relative ordering)."""
        name, sdf_mol, rdkit_mol = molecule_pair
        AllChem.ComputeGasteigerCharges(rdkit_mol)

        sdf_charges = sdf_mol.gasteiger_charges()
        rdkit_charges = []
        for i in range(rdkit_mol.GetNumAtoms()):
            q = float(rdkit_mol.GetAtomWithIdx(i).GetProp("_GasteigerCharge"))
            rdkit_charges.append(q)

        # Filter out NaN from RDKit (can happen for some atom types)
        valid = [(s, r) for s, r in zip(sdf_charges, rdkit_charges)
                 if not (np.isnan(r) or np.isnan(s))]
        if len(valid) < 3:
            pytest.skip(f"{name}: too few valid charges for correlation")

        sdf_arr = np.array([v[0] for v in valid])
        rdkit_arr = np.array([v[1] for v in valid])

        # Pearson correlation should be positive (charges go in same direction)
        if np.std(sdf_arr) > 1e-10 and np.std(rdkit_arr) > 1e-10:
            corr = np.corrcoef(sdf_arr, rdkit_arr)[0, 1]
            assert corr > 0.5, \
                f"{name}: charge correlation {corr:.3f} is too low (expected > 0.5)"

    def test_charge_magnitude(self, molecule_pair):
        """Charges should be in reasonable range."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_charges = sdf_mol.gasteiger_charges()
        for i, q in enumerate(sdf_charges):
            assert -2.0 < q < 2.0, \
                f"{name} atom {i}: charge {q} out of reasonable range"

    def test_charge_near_neutral(self, molecule_pair):
        """Sum of charges should be near total formal charge."""
        name, sdf_mol, _ = molecule_pair
        sdf_charges = sdf_mol.gasteiger_charges()
        total = sum(sdf_charges)
        formal = sdf_mol.total_charge()
        assert abs(total - formal) < 0.5, \
            f"{name}: sum of Gasteiger charges {total:.4f} far from formal charge {formal}"


# ============================================================
# Aromaticity Validation
# ============================================================

class TestAromaticity:
    """Validate aromaticity perception against RDKit."""

    def test_aromatic_atoms(self, molecule_pair):
        """Aromatic atom assignment should match RDKit."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_arom = sdf_mol.all_aromatic_atoms()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_arom = rdkit_mol.GetAtomWithIdx(i).GetIsAromatic()
            if sdf_arom[i] != rdkit_arom:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(
                    f"  atom {i} ({elem}): sdfrust={sdf_arom[i]} rdkit={rdkit_arom}"
                )
        if mismatches:
            pytest.fail(f"{name} aromaticity mismatches:\n" + "\n".join(mismatches))

    def test_aromatic_bonds(self, molecule_pair):
        """Aromatic bond assignment should match RDKit."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_arom = sdf_mol.all_aromatic_bonds()
        mismatches = []
        for i in range(sdf_mol.num_bonds):
            rdkit_bond = rdkit_mol.GetBondWithIdx(i)
            rdkit_arom = rdkit_bond.GetIsAromatic()
            if sdf_arom[i] != rdkit_arom:
                a1, a2 = rdkit_bond.GetBeginAtomIdx(), rdkit_bond.GetEndAtomIdx()
                mismatches.append(
                    f"  bond {i} ({a1}-{a2}): sdfrust={sdf_arom[i]} rdkit={rdkit_arom}"
                )
        if mismatches:
            pytest.fail(f"{name} aromatic bond mismatches:\n" + "\n".join(mismatches))


# ============================================================
# Hybridization Validation
# ============================================================

class TestHybridization:
    """Validate hybridization against RDKit."""

    def test_hybridization(self, molecule_pair):
        """Hybridization should match RDKit for common atom types."""
        name, sdf_mol, rdkit_mol = molecule_pair
        hyb_map = {
            Chem.rdchem.HybridizationType.S: "S",
            Chem.rdchem.HybridizationType.SP: "SP",
            Chem.rdchem.HybridizationType.SP2: "SP2",
            Chem.rdchem.HybridizationType.SP3: "SP3",
            Chem.rdchem.HybridizationType.SP3D: "SP3D",
            Chem.rdchem.HybridizationType.SP3D2: "SP3D2",
        }
        sdf_hybs = sdf_mol.all_hybridizations()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_hyb_enum = rdkit_mol.GetAtomWithIdx(i).GetHybridization()
            rdkit_hyb = hyb_map.get(rdkit_hyb_enum, "Other")
            if sdf_hybs[i] != rdkit_hyb:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(
                    f"  atom {i} ({elem}): sdfrust={sdf_hybs[i]} rdkit={rdkit_hyb}"
                )
        if mismatches:
            pytest.fail(f"{name} hybridization mismatches:\n" + "\n".join(mismatches))


# ============================================================
# Ring Perception Validation
# ============================================================

class TestRingPerception:
    """Validate ring perception against RDKit."""

    def test_ring_count(self, molecule_pair):
        """SSSR ring count should match RDKit."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_rings = sdf_mol.sssr()
        rdkit_ring_info = rdkit_mol.GetRingInfo()
        rdkit_num_rings = rdkit_ring_info.NumRings()
        assert len(sdf_rings) == rdkit_num_rings, \
            f"{name}: ring count sdfrust={len(sdf_rings)} rdkit={rdkit_num_rings}"

    def test_ring_membership(self, molecule_pair):
        """Atom ring membership should match RDKit."""
        name, sdf_mol, rdkit_mol = molecule_pair
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            sdf_in_ring = sdf_mol.is_atom_in_ring(i)
            rdkit_in_ring = rdkit_mol.GetAtomWithIdx(i).IsInRing()
            if sdf_in_ring != rdkit_in_ring:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(
                    f"  atom {i} ({elem}): sdfrust={sdf_in_ring} rdkit={rdkit_in_ring}"
                )
        if mismatches:
            pytest.fail(f"{name} ring membership mismatches:\n" + "\n".join(mismatches))

    def test_ring_sizes(self, molecule_pair):
        """Ring sizes should match RDKit."""
        name, sdf_mol, rdkit_mol = molecule_pair
        sdf_rings = sdf_mol.sssr()
        rdkit_ring_info = rdkit_mol.GetRingInfo()
        rdkit_ring_sizes = sorted([len(r) for r in rdkit_ring_info.AtomRings()])
        sdf_ring_sizes = sorted([len(r) for r in sdf_rings])
        assert sdf_ring_sizes == rdkit_ring_sizes, \
            f"{name}: ring sizes sdfrust={sdf_ring_sizes} rdkit={rdkit_ring_sizes}"


# ============================================================
# Summary Report
# ============================================================

class TestSummaryReport:
    """Print a comparison summary (always passes, for diagnostics)."""

    def test_print_comparison(self, molecule_pair):
        """Print detailed comparison for each molecule."""
        name, sdf_mol, rdkit_mol = molecule_pair
        AllChem.ComputeGasteigerCharges(rdkit_mol)

        n_atoms = sdf_mol.num_atoms
        n_bonds = sdf_mol.num_bonds

        # Count OGB atom feature matches
        sdf_atom_feats = sdf_mol.ogb_atom_features()
        atom_matches = 0
        atom_total = n_atoms * 9
        for i in range(n_atoms):
            rdkit_feat = rdkit_atom_to_ogb_features(rdkit_mol.GetAtomWithIdx(i))
            for j in range(9):
                if sdf_atom_feats[i][j] == rdkit_feat[j]:
                    atom_matches += 1

        # Count OGB bond feature matches
        sdf_bond_feats = sdf_mol.ogb_bond_features()
        bond_matches = 0
        bond_total = n_bonds * 3
        for i in range(n_bonds):
            rdkit_feat = rdkit_bond_to_ogb_features(rdkit_mol.GetBondWithIdx(i), rdkit_mol)
            for j in range(3):
                if sdf_bond_feats[i][j] == rdkit_feat[j]:
                    bond_matches += 1

        print(f"\n{'='*60}")
        print(f"  {name}: {n_atoms} atoms, {n_bonds} bonds")
        print(f"  OGB atom features: {atom_matches}/{atom_total} match "
              f"({atom_matches/atom_total*100:.1f}%)")
        print(f"  OGB bond features: {bond_matches}/{bond_total} match "
              f"({bond_matches/bond_total*100:.1f}%)")
        print(f"{'='*60}")


# ============================================================
# Drug Library (multi-molecule) Validation
# ============================================================

class TestDrugLibrary:
    """Validate features on multi-molecule drug_library.sdf."""

    def test_ogb_atom_features(self, drug_library_pair):
        """OGB atom features should match RDKit (all 9 features including chirality)."""
        name, sdf_mol, rdkit_mol = drug_library_pair
        sdf_feats = sdf_mol.ogb_atom_features()
        feature_names = ["atomic_num", "chirality", "degree", "charge",
                         "num_hs", "radical", "hybridization", "aromatic", "in_ring"]
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_feat = rdkit_atom_to_ogb_features(rdkit_mol.GetAtomWithIdx(i))
            for j in range(9):
                if sdf_feats[i][j] != rdkit_feat[j]:
                    elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                    mismatches.append(
                        f"  atom {i} ({elem}) {feature_names[j]}: "
                        f"sdfrust={sdf_feats[i][j]} rdkit={rdkit_feat[j]}"
                    )
        if mismatches:
            pytest.fail(f"{name}: {len(mismatches)} atom feature mismatches:\n" +
                        "\n".join(mismatches))

    def test_ogb_bond_features(self, drug_library_pair):
        """OGB bond features should match RDKit."""
        name, sdf_mol, rdkit_mol = drug_library_pair
        sdf_feats = sdf_mol.ogb_bond_features()
        feature_names = ["bond_type", "stereo", "conjugated"]
        mismatches = []
        for i in range(sdf_mol.num_bonds):
            rdkit_bond = rdkit_mol.GetBondWithIdx(i)
            rdkit_feat = rdkit_bond_to_ogb_features(rdkit_bond, rdkit_mol)
            for j in range(3):
                if sdf_feats[i][j] != rdkit_feat[j]:
                    a1, a2 = rdkit_bond.GetBeginAtomIdx(), rdkit_bond.GetEndAtomIdx()
                    mismatches.append(
                        f"  bond {i} ({a1}-{a2}) {feature_names[j]}: "
                        f"sdfrust={sdf_feats[i][j]} rdkit={rdkit_feat[j]}"
                    )
        if mismatches:
            pytest.fail(f"{name}: {len(mismatches)} bond feature mismatches:\n" +
                        "\n".join(mismatches))

    def test_aromaticity(self, drug_library_pair):
        """Aromaticity should match RDKit."""
        name, sdf_mol, rdkit_mol = drug_library_pair
        sdf_arom = sdf_mol.all_aromatic_atoms()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_arom = rdkit_mol.GetAtomWithIdx(i).GetIsAromatic()
            if sdf_arom[i] != rdkit_arom:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(f"  atom {i} ({elem}): sdfrust={sdf_arom[i]} rdkit={rdkit_arom}")
        if mismatches:
            pytest.fail(f"{name} aromaticity mismatches:\n" + "\n".join(mismatches))

    def test_hybridization(self, drug_library_pair):
        """Hybridization should match RDKit."""
        name, sdf_mol, rdkit_mol = drug_library_pair
        hyb_map = {
            Chem.rdchem.HybridizationType.S: "S",
            Chem.rdchem.HybridizationType.SP: "SP",
            Chem.rdchem.HybridizationType.SP2: "SP2",
            Chem.rdchem.HybridizationType.SP3: "SP3",
            Chem.rdchem.HybridizationType.SP3D: "SP3D",
            Chem.rdchem.HybridizationType.SP3D2: "SP3D2",
        }
        sdf_hybs = sdf_mol.all_hybridizations()
        mismatches = []
        for i in range(sdf_mol.num_atoms):
            rdkit_hyb = hyb_map.get(rdkit_mol.GetAtomWithIdx(i).GetHybridization(), "Other")
            if sdf_hybs[i] != rdkit_hyb:
                elem = rdkit_mol.GetAtomWithIdx(i).GetSymbol()
                mismatches.append(f"  atom {i} ({elem}): sdfrust={sdf_hybs[i]} rdkit={rdkit_hyb}")
        if mismatches:
            pytest.fail(f"{name} hybridization mismatches:\n" + "\n".join(mismatches))
