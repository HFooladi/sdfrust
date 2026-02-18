#!/usr/bin/env python3
"""
ML Feature Computation Examples for sdfrust Python Bindings

This script demonstrates the molecular ML features:
- OGB-compatible GNN atom and bond featurization
- ECFP/Morgan fingerprints and Tanimoto similarity
- Gasteiger partial charges
- Ring perception (SSSR), aromaticity, hybridization, conjugation
- Valence and hydrogen counts
- NumPy array outputs for direct use with PyTorch Geometric

Run this script from the sdfrust-python directory after building:
    cd sdfrust-python
    maturin develop --features numpy
    python examples/ml_features.py
"""

import os

import sdfrust

# Path to test data
TEST_DATA = os.path.join(os.path.dirname(__file__), "..", "..", "tests", "test_data")


def make_benzene():
    """Create benzene molecule with aromatic bonds."""
    mol = sdfrust.Molecule("benzene")
    for i in range(6):
        mol.add_atom(sdfrust.Atom(i, "C", 0.0, 0.0, 0.0))
    for i in range(6):
        mol.add_bond(sdfrust.Bond(i, (i + 1) % 6, sdfrust.BondOrder.aromatic()))
    return mol


def make_ethanol():
    """Create ethanol molecule (C2H5OH) with coordinates."""
    mol = sdfrust.Molecule("ethanol")
    mol.add_atom(sdfrust.Atom(0, "C", -0.001, 1.086, 0.008))
    mol.add_atom(sdfrust.Atom(1, "C", 0.002, -0.422, 0.002))
    mol.add_atom(sdfrust.Atom(2, "O", 1.210, -0.907, -0.003))
    mol.add_atom(sdfrust.Atom(3, "H", 1.020, 1.465, 0.001))
    mol.add_atom(sdfrust.Atom(4, "H", -0.536, 1.449, -0.873))
    mol.add_atom(sdfrust.Atom(5, "H", -0.523, 1.437, 0.902))
    mol.add_atom(sdfrust.Atom(6, "H", -0.524, -0.785, 0.891))
    mol.add_atom(sdfrust.Atom(7, "H", -0.510, -0.793, -0.884))
    mol.add_atom(sdfrust.Atom(8, "H", 1.192, -1.870, -0.003))
    mol.add_bond(sdfrust.Bond(0, 1, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 2, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 3, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 4, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(0, 5, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 6, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(1, 7, sdfrust.BondOrder.single()))
    mol.add_bond(sdfrust.Bond(2, 8, sdfrust.BondOrder.single()))
    return mol


def example_ring_perception():
    """Demonstrate SSSR ring perception."""
    print("=" * 60)
    print("SSSR Ring Perception")
    print("=" * 60)

    mol = make_benzene()
    rings = mol.sssr()
    print(f"Benzene has {len(rings)} ring(s)")
    for i, ring in enumerate(rings):
        print(f"  Ring {i}: atoms {ring} (size {len(ring)})")

    # Ring sizes for atoms
    for i in range(min(3, mol.num_atoms)):
        sizes = mol.ring_sizes(i)
        smallest = mol.smallest_ring_size(i)
        print(f"  Atom {i}: ring sizes = {sizes}, smallest = {smallest}")
    print()


def example_aromaticity():
    """Demonstrate aromaticity detection."""
    print("=" * 60)
    print("Aromaticity Detection")
    print("=" * 60)

    mol = make_benzene()
    print("Benzene aromaticity:")
    aromatic_atoms = mol.all_aromatic_atoms()
    aromatic_bonds = mol.all_aromatic_bonds()
    print(f"  Aromatic atoms: {aromatic_atoms}")
    print(f"  Aromatic bonds: {aromatic_bonds}")

    # Per-atom check
    for i in range(mol.num_atoms):
        print(f"  Atom {i}: is_aromatic = {mol.is_aromatic_atom(i)}")
    print()


def example_hybridization():
    """Demonstrate hybridization inference."""
    print("=" * 60)
    print("Hybridization")
    print("=" * 60)

    mol = make_ethanol()
    hybs = mol.all_hybridizations()
    for i, hyb in enumerate(hybs):
        atom = mol.atoms[i]
        print(f"  Atom {i} ({atom.element}): {hyb}")
    print()


def example_conjugation():
    """Demonstrate conjugation detection."""
    print("=" * 60)
    print("Conjugation Detection")
    print("=" * 60)

    mol = make_benzene()
    conj = mol.all_conjugated_bonds()
    print("Benzene bond conjugation:")
    for i, is_conj in enumerate(conj):
        bond = mol.bonds[i]
        print(f"  Bond {i} ({bond.atom1}-{bond.atom2}): conjugated = {is_conj}")
    print()


def example_valence():
    """Demonstrate valence and hydrogen counts."""
    print("=" * 60)
    print("Valence & Hydrogen Counts")
    print("=" * 60)

    mol = make_ethanol()
    for i in range(min(3, mol.num_atoms)):
        atom = mol.atoms[i]
        degree = mol.atom_degree(i)
        total_h = mol.total_hydrogen_count(i)
        implicit_h = mol.implicit_hydrogen_count(i)
        print(
            f"  Atom {i} ({atom.element}): degree={degree}, "
            f"total_H={total_h}, implicit_H={implicit_h}"
        )
    print()


def example_ogb_features():
    """Demonstrate OGB-compatible GNN featurization."""
    print("=" * 60)
    print("OGB GNN Features")
    print("=" * 60)

    mol = make_benzene()

    # Atom features [N, 9]
    atom_feats = mol.ogb_atom_features()
    print(f"Atom features: {len(atom_feats)} atoms x {len(atom_feats[0])} features")
    print("Feature names: [atomic_num, chirality, degree, charge+5, "
          "num_hs, radical, hybridization, is_aromatic, is_in_ring]")
    for i, feat in enumerate(atom_feats[:3]):
        print(f"  Atom {i}: {feat}")
    if len(atom_feats) > 3:
        print(f"  ... ({len(atom_feats) - 3} more)")

    # Bond features [E, 3]
    bond_feats = mol.ogb_bond_features()
    print(f"\nBond features: {len(bond_feats)} bonds x {len(bond_feats[0])} features")
    print("Feature names: [bond_type, stereo, is_conjugated]")
    for i, feat in enumerate(bond_feats[:3]):
        print(f"  Bond {i}: {feat}")

    # Full graph features (dict)
    graph = mol.ogb_graph_features()
    print(f"\nFull graph:")
    print(f"  Atoms: {graph['num_atoms']}")
    print(f"  Directed edges: {graph['num_bonds']}")
    print(f"  Edge src (first 6): {graph['edge_src'][:6]}")
    print(f"  Edge dst (first 6): {graph['edge_dst'][:6]}")
    print()


def example_ecfp_fingerprints():
    """Demonstrate ECFP/Morgan fingerprints."""
    print("=" * 60)
    print("ECFP/Morgan Fingerprints")
    print("=" * 60)

    benzene = make_benzene()
    ethanol = make_ethanol()

    # Compute ECFP4 (radius=2)
    fp_benzene = benzene.ecfp(radius=2, n_bits=2048)
    fp_ethanol = ethanol.ecfp(radius=2, n_bits=2048)
    print(f"Benzene ECFP4: {sum(fp_benzene)} bits set / {len(fp_benzene)}")
    print(f"Ethanol ECFP4: {sum(fp_ethanol)} bits set / {len(fp_ethanol)}")

    # On-bit indices
    on_bits = benzene.ecfp_on_bits(radius=2, n_bits=2048)
    print(f"Benzene on-bits: {on_bits[:10]}...")

    # Tanimoto similarity
    sim_self = benzene.tanimoto_similarity(benzene)
    sim_other = benzene.tanimoto_similarity(ethanol)
    print(f"\nTanimoto(benzene, benzene): {sim_self:.4f}")
    print(f"Tanimoto(benzene, ethanol): {sim_other:.4f}")

    # Count fingerprint
    counts = benzene.ecfp_counts(radius=2)
    print(f"\nBenzene count FP: {len(counts)} unique features")
    for hash_val, count in list(counts.items())[:5]:
        print(f"  Feature {hash_val}: count={count}")
    print()


def example_gasteiger_charges():
    """Demonstrate Gasteiger partial charges."""
    print("=" * 60)
    print("Gasteiger Partial Charges")
    print("=" * 60)

    mol = make_ethanol()
    charges = mol.gasteiger_charges()
    for i, charge in enumerate(charges):
        atom = mol.atoms[i]
        print(f"  Atom {i} ({atom.element}): {charge:+.4f}")

    total = sum(charges)
    print(f"  Sum of charges: {total:+.6f}")
    print()


def example_numpy_ml_features():
    """Demonstrate NumPy array outputs for ML pipelines."""
    print("=" * 60)
    print("NumPy ML Feature Arrays")
    print("=" * 60)

    try:
        import numpy as np
    except ImportError:
        print("NumPy not available, skipping this example")
        print()
        return

    mol = make_benzene()

    # OGB atom features as NumPy array
    atom_feats = mol.get_ogb_atom_features_array()
    print(f"OGB atom features: shape={atom_feats.shape}, dtype={atom_feats.dtype}")
    print(f"  First row: {atom_feats[0]}")

    # OGB bond features as NumPy array
    bond_feats = mol.get_ogb_bond_features_array()
    print(f"OGB bond features: shape={bond_feats.shape}, dtype={bond_feats.dtype}")

    # Gasteiger charges as NumPy array
    charges = mol.get_gasteiger_charges_array()
    print(f"Gasteiger charges: shape={charges.shape}, dtype={charges.dtype}")
    print(f"  Values: {charges}")

    # ECFP as NumPy array
    fp = mol.get_ecfp_array(radius=2, n_bits=2048)
    print(f"ECFP fingerprint: shape={fp.shape}, dtype={fp.dtype}")
    print(f"  Bits set: {np.sum(fp)}")

    # Coordinates
    coords = mol.get_coords_array()
    print(f"Coordinates: shape={coords.shape}")

    print("\n--- Ready for PyTorch Geometric ---")
    print("  atom_features = torch.tensor(mol.get_ogb_atom_features_array())")
    print("  bond_features = torch.tensor(mol.get_ogb_bond_features_array())")
    print("  graph = mol.ogb_graph_features()")
    print("  edge_index = torch.tensor([graph['edge_src'], graph['edge_dst']])")
    print()


def example_real_molecule():
    """Demonstrate ML features on a real molecule from file."""
    print("=" * 60)
    print("Real Molecule ML Features")
    print("=" * 60)

    aspirin_path = os.path.join(TEST_DATA, "aspirin.sdf")
    if not os.path.exists(aspirin_path):
        print("Aspirin test file not found, skipping")
        print()
        return

    mol = sdfrust.parse_sdf_file(aspirin_path)
    print(f"Molecule: {mol.name}")
    print(f"Formula: {mol.formula()}")
    print(f"Atoms: {mol.num_atoms}, Bonds: {mol.num_bonds}")

    # Rings
    rings = mol.sssr()
    print(f"\nSSSR rings: {len(rings)}")
    for i, ring in enumerate(rings):
        print(f"  Ring {i}: size {len(ring)}")

    # Aromaticity
    aromatic = mol.all_aromatic_atoms()
    n_aromatic = sum(aromatic)
    print(f"\nAromatic atoms: {n_aromatic} / {mol.num_atoms}")

    # Hybridization distribution
    hybs = mol.all_hybridizations()
    hyb_counts = {}
    for h in hybs:
        hyb_counts[h] = hyb_counts.get(h, 0) + 1
    print(f"Hybridization distribution: {hyb_counts}")

    # OGB features
    atom_feats = mol.ogb_atom_features()
    bond_feats = mol.ogb_bond_features()
    print(f"\nOGB features: [{len(atom_feats)}, 9] atoms, [{len(bond_feats)}, 3] bonds")

    # ECFP
    fp = mol.ecfp(radius=2, n_bits=2048)
    print(f"ECFP4: {sum(fp)} / {len(fp)} bits set ({sum(fp)/len(fp)*100:.1f}% density)")

    # Gasteiger charges
    charges = mol.gasteiger_charges()
    min_q = min(charges)
    max_q = max(charges)
    print(f"Gasteiger charges: min={min_q:+.4f}, max={max_q:+.4f}")
    print()


def main():
    """Run all ML feature examples."""
    print("\n" + "=" * 60)
    print("  sdfrust ML Feature Examples")
    print("=" * 60 + "\n")

    example_ring_perception()
    example_aromaticity()
    example_hybridization()
    example_conjugation()
    example_valence()
    example_ogb_features()
    example_ecfp_fingerprints()
    example_gasteiger_charges()
    example_numpy_ml_features()
    example_real_molecule()

    print("=" * 60)
    print("All ML feature examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
