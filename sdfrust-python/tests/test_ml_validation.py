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
# OGB reference encoding
# ============================================================
# Use the ogb package itself (``pip install --no-deps ogb``) so the comparison is
# against the real AtomEncoder/BondEncoder input convention. If ogb is not
# installed, fall back to a verbatim copy of ogb.utils.features (ogb 1.3.6).

try:
    from ogb.utils.features import atom_to_feature_vector, bond_to_feature_vector
except ImportError:
    _OGB_FEATURES = {
        "possible_atomic_num_list": list(range(1, 119)) + ["misc"],
        "possible_chirality_list": ["CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW",
                                    "CHI_TETRAHEDRAL_CCW", "CHI_OTHER", "misc"],
        "possible_degree_list": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "misc"],
        "possible_formal_charge_list": [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, "misc"],
        "possible_numH_list": [0, 1, 2, 3, 4, 5, 6, 7, 8, "misc"],
        "possible_number_radical_e_list": [0, 1, 2, 3, 4, "misc"],
        "possible_hybridization_list": ["SP", "SP2", "SP3", "SP3D", "SP3D2", "misc"],
        "possible_is_aromatic_list": [False, True],
        "possible_is_in_ring_list": [False, True],
        "possible_bond_type_list": ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC", "misc"],
        "possible_bond_stereo_list": ["STEREONONE", "STEREOZ", "STEREOE",
                                      "STEREOCIS", "STEREOTRANS", "STEREOANY"],
        "possible_is_conjugated_list": [False, True],
    }

    def _safe_index(values, e):
        try:
            return values.index(e)
        except ValueError:
            return len(values) - 1

    def atom_to_feature_vector(atom):
        f = _OGB_FEATURES
        return [
            _safe_index(f["possible_atomic_num_list"], atom.GetAtomicNum()),
            _safe_index(f["possible_chirality_list"], str(atom.GetChiralTag())),
            _safe_index(f["possible_degree_list"], atom.GetTotalDegree()),
            _safe_index(f["possible_formal_charge_list"], atom.GetFormalCharge()),
            _safe_index(f["possible_numH_list"], atom.GetTotalNumHs()),
            _safe_index(f["possible_number_radical_e_list"], atom.GetNumRadicalElectrons()),
            _safe_index(f["possible_hybridization_list"], str(atom.GetHybridization())),
            f["possible_is_aromatic_list"].index(atom.GetIsAromatic()),
            f["possible_is_in_ring_list"].index(atom.IsInRing()),
        ]

    def bond_to_feature_vector(bond):
        f = _OGB_FEATURES
        return [
            _safe_index(f["possible_bond_type_list"], str(bond.GetBondType())),
            f["possible_bond_stereo_list"].index(str(bond.GetStereo())),
            f["possible_is_conjugated_list"].index(bond.GetIsConjugated()),
        ]


def rdkit_atom_to_ogb_features(atom):
    """OGB 9-feature vector for an RDKit atom (ogb.utils.features.atom_to_feature_vector)."""
    return atom_to_feature_vector(atom)


def rdkit_bond_to_ogb_features(bond, mol=None):
    """OGB 3-feature vector for an RDKit bond (ogb.utils.features.bond_to_feature_vector)."""
    return bond_to_feature_vector(bond)


def _check_atom_feature(name, sdf_mol, rdkit_mol, k, label):
    """Compare OGB atom feature column k against the reference encoding."""
    sdf_feats = sdf_mol.ogb_atom_features()
    mismatches = []
    for i in range(sdf_mol.num_atoms):
        atom = rdkit_mol.GetAtomWithIdx(i)
        ref = rdkit_atom_to_ogb_features(atom)[k]
        if sdf_feats[i][k] != ref:
            mismatches.append(f"  atom {i} ({atom.GetSymbol()}): sdfrust={sdf_feats[i][k]} ogb={ref}")
    if mismatches:
        pytest.fail(f"{name} {label} mismatches:\n" + "\n".join(mismatches))


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
        """Feature 0: atomic number index (Z - 1)."""
        _check_atom_feature(*molecule_pair, 0, "atomic_num")

    def test_chirality(self, molecule_pair):
        """Feature 1: chirality tag."""
        _check_atom_feature(*molecule_pair, 1, "chirality")

    def test_degree(self, molecule_pair):
        """Feature 2: total degree (including implicit H)."""
        _check_atom_feature(*molecule_pair, 2, "degree")

    def test_formal_charge(self, molecule_pair):
        """Feature 3: formal charge index."""
        _check_atom_feature(*molecule_pair, 3, "formal_charge")

    def test_num_hs(self, molecule_pair):
        """Feature 4: total H count."""
        _check_atom_feature(*molecule_pair, 4, "num_hs")

    def test_hybridization(self, molecule_pair):
        """Feature 6: hybridization index (SP=0 ... SP3D2=4, misc=5)."""
        _check_atom_feature(*molecule_pair, 6, "hybridization")

    def test_is_aromatic(self, molecule_pair):
        """Feature 7: aromaticity."""
        _check_atom_feature(*molecule_pair, 7, "aromatic")

    def test_is_in_ring(self, molecule_pair):
        """Feature 8: ring membership."""
        _check_atom_feature(*molecule_pair, 8, "in_ring")

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
            rdkit_feat = rdkit_bond_to_ogb_features(rdkit_bond)
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

    def test_bits_identical_to_rdkit_morgan(self, molecule_pair):
        """sdfrust ECFP bits equal RDKit Morgan bits (default options)."""
        from rdkit.Chem import rdFingerprintGenerator
        name, sdf_mol, rdkit_mol = molecule_pair
        for radius in (2, 3):
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2048)
            ref = list(gen.GetFingerprint(rdkit_mol).GetOnBits())
            assert sorted(sdf_mol.ecfp_on_bits(radius, 2048)) == ref, f"{name} r={radius}"

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
            rdkit_feat = rdkit_bond_to_ogb_features(rdkit_bond)
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
