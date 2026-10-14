#!/usr/bin/env python
"""Differential comparison of sdfrust against RDKit, OpenBabel and the OGB package.

Runs every sdfrust feature that has an equivalent in an established toolkit on the
same input files and reports per-feature agreement plus the disagreeing rows.

Usage:
    python benchmarks/compare_tools.py [--out DIR] [--label NAME] [--pdbbind N]

Requires: rdkit, ogb (``pip install --no-deps ogb``), optionally the ``obabel`` CLI.
Writes ``<out>/<label>.md`` and ``<out>/<label>.json``.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import sdfrust
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

try:
    from ogb.utils.features import atom_to_feature_vector, bond_to_feature_vector
except ImportError:  # pragma: no cover
    atom_to_feature_vector = bond_to_feature_vector = None

ROOT = Path(__file__).resolve().parents[2]
OBABEL = shutil.which("obabel")

# Diverse molecules covering the chemistry edge cases found in the audit.
SMILES = {
    "naphthalene": "c1ccc2ccccc2c1",
    "anthracene": "c1ccc2cc3ccccc3cc2c1",
    "phenanthrene": "c1ccc2c(c1)ccc1ccccc12",
    "indole": "c1ccc2[nH]ccc2c1",
    "azulene": "c1ccc2cccc2cc1",
    "quinoline": "c1ccc2ncccc2c1",
    "purine": "c1ncc2[nH]cnc2n1",
    "pyrrole": "c1cc[nH]c1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "imidazole": "c1c[nH]cn1",
    "2-pyridone": "O=c1cccc[nH]1",
    "tropylium": "c1cc[cH+]ccc1",
    "cyclopentadienide": "[cH-]1cccc1",
    "pyridinium": "c1cc[nH+]cc1",
    "biphenyl": "c1ccc(-c2ccccc2)cc1",
    "dimethylsulfone": "CS(C)(=O)=O",
    "dmso": "CS(C)=O",
    "phosphate": "COP(=O)(OC)OC",
    "allene": "C=C=C",
    "co2": "O=C=O",
    "azide": "CN=[N+]=[N-]",
    "acetonitrile": "CC#N",
    "benzonitrile": "N#Cc1ccccc1",
    "nitrobenzene": "O=[N+]([O-])c1ccccc1",
    "acetamide": "CC(N)=O",
    "nmethylacetamide": "CNC(C)=O",
    "aniline": "Nc1ccccc1",
    "anilinium": "[NH3+]c1ccccc1",
    "phenol": "Oc1ccccc1",
    "methyl_acetate": "COC(C)=O",
    "acetate": "CC(=O)[O-]",
    "glycine_zw": "[NH3+]CC(=O)[O-]",
    "tetramethylammonium": "C[N+](C)(C)C",
    "butadiene": "C=CC=C",
    "enyne": "C=CC#C",
    "cubane": "C12C3C4C1C1C2C3C41",
    "adamantane": "C1C2CC3CC1CC(C2)C3",
    "norbornane": "C1CC2CCC1C2",
    "spiro": "C1CCC2(C1)CCCC2",
    "l_alanine": "C[C@H](N)C(=O)O",
    "d_alanine": "C[C@@H](N)C(=O)O",
    "r_bromochlorofluoromethane": "F[C@H](Cl)Br",
    "e_2_butene": "C/C=C/C",
    "z_2_butene": "C/C=C\\C",
    "c13_methane": "[13CH4]",
    "deuterated_benzene": "[2H]c1ccccc1",
    "uranium235_f6": "F[235U](F)(F)(F)(F)F",
    "methyl_radical": "[CH3]",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "sildenafil": "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
    "atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
    "paclitaxel_core": "CC1=C2[C@@H](OC(C)=O)C(=O)[C@]3(C)[C@@H](O)C[C@H]4OC[C@@]4(OC(C)=O)[C@H]3[C@H](OC(=O)c3ccccc3)[C@](O)(C[C@@H]1O)C2(C)C",
    "cisplatin": "N.N.Cl[Pt]Cl",
    "ferrocene_like": "[Fe+2]",
    "selenophene": "c1cc[se]c1",
    "boronic_acid": "OB(O)c1ccccc1",
    "trifluoromethylbenzene": "FC(F)(F)c1ccccc1",
    "iodobenzene": "Ic1ccccc1",
}


# ----------------------------------------------------------------------------
# Corpus
# ----------------------------------------------------------------------------


def build_smiles_corpus(gen: Path, smiles_file: Path, limit: int) -> list[Path]:
    """Convert a SMILES file into two multi-record Kekulé SDFs (heavy-atom only, explicit H)."""
    noh = gen / f"{smiles_file.stem}__noh.sdf"
    withh = gen / f"{smiles_file.stem}__h.sdf"
    w_noh, w_h = Chem.SDWriter(str(noh)), Chem.SDWriter(str(withh))
    w_noh.SetKekulize(True)
    w_h.SetKekulize(True)
    n = 0
    for line in smiles_file.read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].upper() == "SMILES":
            continue
        mol = Chem.MolFromSmiles(parts[0])
        if mol is None or mol.GetNumAtoms() == 0:
            continue
        mol.SetProp("_Name", parts[1] if len(parts) > 1 else f"mol{n}")
        AllChem.Compute2DCoords(mol)
        w_noh.write(mol)
        molh = Chem.AddHs(mol, addCoords=True)
        w_h.write(molh)
        n += 1
        if n >= limit:
            break
    w_noh.close()
    w_h.close()
    return [noh, withh]


def build_corpus(workdir: Path, pdbbind_n: int, smiles: list[Path] = (), smiles_limit: int = 2000,
                 extra: list[Path] = ()) -> list[Path]:
    """Return SDF files to compare: repo test data + generated edge cases."""
    files = sorted((ROOT / "tests/test_data").glob("*.sdf"))
    files += sorted((ROOT / "sdfrust-python/examples/data").glob("*.sdf"))

    gen = workdir / "generated"
    gen.mkdir(parents=True, exist_ok=True)
    for name, smi in SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        molh = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(molh, randomSeed=42) != 0:
            AllChem.Compute2DCoords(molh)
        Chem.AssignStereochemistryFrom3D(molh)
        molh.SetProp("_Name", name)
        # Kekulé V2000 with explicit H (the PubChem/ChEMBL 3D style)
        Chem.Kekulize(molh, clearAromaticFlags=True)
        (gen / f"{name}__h.sdf").write_text(Chem.MolToMolBlock(molh, kekulize=True) + "$$$$\n")
        # Heavy-atom-only variant (the typical 2D database style)
        heavy = Chem.RemoveHs(Chem.MolFromSmiles(smi))
        AllChem.Compute2DCoords(heavy)
        heavy.SetProp("_Name", name)
        (gen / f"{name}__noh.sdf").write_text(Chem.MolToMolBlock(heavy, kekulize=True) + "$$$$\n")
    files += sorted(gen.glob("*.sdf"))
    for smi in smiles:
        files += build_smiles_corpus(gen, Path(smi), smiles_limit)
    files += [Path(p) for p in extra]

    pdb_dir = os.environ.get("PDBBIND_2024_DIR")
    if pdbbind_n and pdb_dir:
        ligs = sorted(Path(pdb_dir).rglob("*_ligand.sdf"))[:pdbbind_n]
        files += ligs
    return files


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


class Tally:
    def __init__(self):
        self.ok = defaultdict(int)
        self.total = defaultdict(int)
        self.diffs = defaultdict(list)
        self.numeric = defaultdict(list)
        self.worst = {}

    def check(self, feature, ours, theirs, ctx, max_diffs=40):
        self.total[feature] += 1
        if ours == theirs:
            self.ok[feature] += 1
        elif len(self.diffs[feature]) < max_diffs:
            self.diffs[feature].append({"ctx": ctx, "sdfrust": ours, "reference": theirs})

    def close(self, feature, ours, theirs, ctx, tol):
        self.numeric[feature].append(abs(ours - theirs))
        worst = self.worst.setdefault(feature, [])
        worst.append((abs(ours - theirs), ctx, round(ours, 4), round(theirs, 4)))
        worst.sort(reverse=True)
        del worst[5:]
        self.check(feature, True, abs(ours - theirs) <= tol, ctx)
        if abs(ours - theirs) > tol and self.diffs[feature] and self.diffs[feature][-1]["ctx"] == ctx:
            self.diffs[feature][-1].update({"sdfrust": round(ours, 4), "reference": round(theirs, 4)})

    def crash(self, feature, ctx, exc):
        self.total[feature] += 1
        if len(self.diffs[feature]) < 40:
            self.diffs[feature].append({"ctx": ctx, "sdfrust": f"ERROR {type(exc).__name__}: {exc}"[:200], "reference": "ok"})


HYB = {"S": "S", "SP": "SP", "SP2": "SP2", "SP3": "SP3", "SP3D": "SP3D", "SP3D2": "SP3D2"}


def rdkit_hyb(atom):
    return HYB.get(str(atom.GetHybridization()), "Other")


def rdkit_chiral(atom):
    return {
        Chem.ChiralType.CHI_UNSPECIFIED: "Unspecified",
        Chem.ChiralType.CHI_TETRAHEDRAL_CW: "CW",
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW: "CCW",
    }.get(atom.GetChiralTag(), "Other")


def obabel_props(path: Path):
    """MW / exact mass / formula / Gasteiger charges from OpenBabel."""
    if not OBABEL:
        return None
    try:
        out = subprocess.run(
            [OBABEL, str(path), "-l", "1", "-otxt", "--append", "MW formula"],
            capture_output=True, text=True, timeout=30,
        ).stdout.split()
        mw, formula = float(out[-2]), out[-1]
        mol2 = subprocess.run(
            [OBABEL, str(path), "-l", "1", "-omol2", "--partialcharge", "gasteiger"],
            capture_output=True, text=True, timeout=30,
        ).stdout
        charges, in_atoms = [], False
        for line in mol2.splitlines():
            if line.startswith("@<TRIPOS>"):
                in_atoms = line.strip() == "@<TRIPOS>ATOM"
                continue
            if in_atoms and line.strip():
                charges.append(float(line.split()[8]))
        return {"mw": mw, "formula": formula, "gasteiger": charges}
    except Exception:
        return None


def tanimoto(a, b):
    inter = sum(1 for x, y in zip(a, b) if x and y)
    union = sum(1 for x, y in zip(a, b) if x or y)
    return inter / union if union else 1.0


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    sy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return cov / (sx * sy) if sx and sy else float("nan")


# ----------------------------------------------------------------------------
# Per-molecule comparison
# ----------------------------------------------------------------------------


def compare_molecule(t: Tally, path: Path, idx: int, ours, rd, ob, raw=None):
    ctx = f"{path.name}#{idx}"
    if rd is None:
        t.check("rdkit_can_parse", True, False, ctx)
        return
    n = rd.GetNumAtoms()

    t.check("parse.atom_count", ours.num_atoms, n, ctx)
    t.check("parse.bond_count", ours.num_bonds, rd.GetNumBonds(), ctx)
    if ours.num_atoms != n:
        return
    ours_atoms = ours.atoms
    t.check("parse.formal_charges", [a.formal_charge for a in ours_atoms],
            [a.GetFormalCharge() for a in rd.GetAtoms()], ctx)
    our_iso = [getattr(a, "isotope", None) or (a.mass_difference or None) for a in ours_atoms]
    t.check("parse.isotopes", our_iso, [a.GetIsotope() or None for a in rd.GetAtoms()], ctx)
    if hasattr(ours_atoms[0], "radical_electrons"):
        t.check("parse.radical_electrons", [a.radical_electrons for a in ours_atoms],
                [a.GetNumRadicalElectrons() for a in rd.GetAtoms()], ctx)

    # Molecular descriptors
    _safe(t, "desc.molecular_weight", ctx, lambda: t.close(
        "desc.molecular_weight", ours.molecular_weight() or 0.0, Descriptors.MolWt(rd), ctx, 0.05))
    _safe(t, "desc.exact_mass", ctx, lambda: t.close(
        "desc.exact_mass", ours.exact_mass() or 0.0, Descriptors.ExactMolWt(rd), ctx, 0.005))
    _safe(t, "desc.formula", ctx, lambda: t.check(
        "desc.formula", ours.formula(), rdMolDescriptors.CalcMolFormula(rd), ctx))
    _safe(t, "desc.ring_count", ctx, lambda: t.check(
        "desc.ring_count", ours.ring_count(), len(Chem.GetSSSR(rd)), ctx))
    _safe(t, "desc.rotatable_bonds", ctx, lambda: t.check(
        "desc.rotatable_bonds", ours.rotatable_bond_count(), rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(rd)), ctx))
    _safe(t, "desc.sssr_sizes", ctx, lambda: t.check(
        "desc.sssr_sizes", sorted(len(r) for r in ours.sssr()), sorted(len(r) for r in Chem.GetSSSR(rd)), ctx))

    if ob:
        t.close("openbabel.molecular_weight", ours.molecular_weight() or 0.0, ob["mw"], ctx, 0.05)
        t.check("openbabel.formula", ours.formula(), ob["formula"], ctx)

    # Per-atom perception
    def per_atom():
        arom = ours.all_aromatic_atoms()
        hyb = ours.all_hybridizations()
        chi = ours.all_chiralities()
        for i, a in enumerate(rd.GetAtoms()):
            actx = f"{ctx} atom {i} {a.GetSymbol()}"
            t.check("atom.aromatic", bool(arom[i]), a.GetIsAromatic(), actx)
            t.check("atom.hybridization", hyb[i], rdkit_hyb(a), actx)
            t.check("atom.implicit_h", ours.implicit_hydrogen_count(i), a.GetTotalNumHs(includeNeighbors=False), actx)
            t.check("atom.total_degree", ours.atom_degree(i) + ours.implicit_hydrogen_count(i), a.GetTotalDegree(), actx)
            t.check("atom.in_ring", bool(ours.is_atom_in_ring(i)), a.IsInRing(), actx)
            t.check("atom.chiral_tag", chi[i], rdkit_chiral(a), actx)

    _safe(t, "atom.perception", ctx, per_atom)

    def per_bond():
        conj = ours.all_conjugated_bonds()
        arom = ours.all_aromatic_bonds()
        for b in rd.GetBonds():
            j = b.GetIdx()
            bctx = f"{ctx} bond {j} {b.GetBeginAtomIdx()}-{b.GetEndAtomIdx()}"
            t.check("bond.conjugated", bool(conj[j]), b.GetIsConjugated(), bctx)
            t.check("bond.aromatic", bool(arom[j]), b.GetIsAromatic(), bctx)
            t.check("bond.in_ring", bool(ours.is_bond_in_ring(j)), b.IsInRing(), bctx)

    _safe(t, "bond.perception", ctx, per_bond)

    # CIP R/S (only if sdfrust exposes it)
    if hasattr(ours, "cip_labels"):
        rd_cip = dict(Chem.FindMolChiralCenters(rd, includeUnassigned=False, useLegacyImplementation=False))
        ours_cip = {i: l for i, l in enumerate(ours.cip_labels()) if l}
        t.check("atom.cip_label", ours_cip, rd_cip, ctx)

    # Gasteiger
    def gasteiger():
        rdc = Chem.Mol(rd)
        AllChem.ComputeGasteigerCharges(rdc)
        ref = [a.GetDoubleProp("_GasteigerCharge") for a in rdc.GetAtoms()]
        if any(math.isnan(x) for x in ref):
            return
        mine = ours.gasteiger_charges()
        for i, (x, y) in enumerate(zip(mine, ref)):
            t.close("charge.gasteiger_vs_rdkit", x, y, f"{ctx} atom {i}", 0.02)
        if ob and len(ob["gasteiger"]) == len(mine):
            for i, (x, y) in enumerate(zip(mine, ob["gasteiger"])):
                t.close("charge.gasteiger_vs_openbabel", x, y, f"{ctx} atom {i}", 0.02)

    _safe(t, "charge.gasteiger_vs_rdkit", ctx, gasteiger)

    # OGB encodings vs the ogb package itself
    if atom_to_feature_vector is not None:
        def ogb():
            ours_a = ours.ogb_atom_features()
            names = ["atomic_num", "chirality", "degree", "formal_charge", "num_h",
                     "radical_e", "hybridization", "aromatic", "in_ring"]
            for i, a in enumerate(rd.GetAtoms()):
                ref = atom_to_feature_vector(a)
                for k, name in enumerate(names):
                    t.check(f"ogb.atom.{name}", list(ours_a[i])[k], ref[k], f"{ctx} atom {i} {a.GetSymbol()}")
            ours_b = ours.ogb_bond_features()
            bnames = ["bond_type", "stereo", "conjugated"]
            for b in rd.GetBonds():
                ref = bond_to_feature_vector(b)
                for k, name in enumerate(bnames):
                    t.check(f"ogb.bond.{name}", list(ours_b[b.GetIdx()])[k], ref[k], f"{ctx} bond {b.GetIdx()}")

        _safe(t, "ogb.atom", ctx, ogb)

    # ECFP vs RDKit Morgan: bits should be identical on the same graph
    def ecfp_bits():
        for radius in (2, 3):
            gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2048)
            t.check(f"ecfp.morgan_bits_r{radius}", sorted(ours.ecfp_on_bits(radius, 2048)),
                    list(gen.GetFingerprint(rd).GetOnBits()), ctx)

    _safe(t, "ecfp.morgan_bits_r2", ctx, ecfp_bits)

    # Round trip: sdfrust writer -> RDKit reader
    def roundtrip():
        block = sdfrust.write_sdf_auto_string(ours)
        back = Chem.MolFromMolBlock(block.split("$$$$")[0], sanitize=False, removeHs=False)
        if back is None:
            t.check("roundtrip.rdkit_reads_our_output", False, True, ctx)
            return
        t.check("roundtrip.atom_count", back.GetNumAtoms(), n, ctx)
        t.check("roundtrip.formal_charges", [a.GetFormalCharge() for a in back.GetAtoms()],
                [a.GetFormalCharge() for a in rd.GetAtoms()], ctx)
        t.check("roundtrip.isotopes", [a.GetIsotope() for a in back.GetAtoms()],
                [a.GetIsotope() for a in rd.GetAtoms()], ctx)
        t.check("roundtrip.radicals", [a.GetNumRadicalElectrons() for a in back.GetAtoms()],
                [a.GetNumRadicalElectrons() for a in rd.GetAtoms()], ctx)
        t.check("roundtrip.bond_orders", [str(b.GetBondType()) for b in back.GetBonds()],
                [str(b.GetBondType()) for b in (raw or rd).GetBonds()], ctx)

    _safe(t, "roundtrip", ctx, roundtrip)


def _safe(t, feature, ctx, fn):
    try:
        fn()
    except BaseException as exc:  # PanicException derives from BaseException
        if isinstance(exc, KeyboardInterrupt):
            raise
        t.crash(feature, ctx, exc)


# ----------------------------------------------------------------------------
# Targeted parser/writer robustness checks (each must hold for a correct library)
# ----------------------------------------------------------------------------


def edge_checks(workdir: Path) -> list[dict]:
    results = []

    def run(name, fn, expect):
        try:
            got = fn()
            ok = got == expect
            results.append({"check": name, "ok": ok, "got": repr(got)[:160], "expected": repr(expect)})
        except BaseException as exc:
            if isinstance(exc, KeyboardInterrupt):
                raise
            results.append({"check": name, "ok": False, "got": f"{type(exc).__name__}: {exc}"[:160], "expected": repr(expect)})

    v3000 = (ROOT / "tests/test_data/v3000_benzene.sdf").read_text()
    v2000 = (ROOT / "tests/test_data/aspirin.sdf").read_text()
    run("V3000 via parse_sdf_string gives atoms (or error), not empty", lambda: sdfrust.parse_sdf_string(v3000).num_atoms, 6)
    run("Mixed V2000+V3000 file parses both records", lambda: [m.num_atoms for m in sdfrust.parse_auto_string_multi(v2000 + v3000)], [21, 6])
    run("Trailing blank line after $$$$", lambda: len(sdfrust.parse_sdf_string_multi(v2000 + "\n\n")), 1)
    bad = "bad\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\nGARBAGE\nM  END\n$$$$\n"
    run("Iterator recovers after malformed record", lambda: _iter_count(workdir, bad + v2000 + v2000), 2)
    nonascii = v2000.replace("\n", "\nµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅµÅ\n", 4).split("\n", 4)
    run("Non-ASCII junk raises an error, not a panic", lambda: _raises_clean(lambda: sdfrust.parse_sdf_string("\n".join(nonascii))), True)
    huge = "big\n\n\n  1  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 999999999999 0 0 0 0\nM  V30 END CTAB\nM  END\n"
    run("Huge V3000 count raises an error, not an abort", lambda: _raises_clean(lambda: sdfrust.parse_sdf_v3000_string(huge)), True)

    mol2_am = """@<TRIPOS>MOLECULE
acetamide
 4 3 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1  0.0 0.0 0.0 C.3  1 LIG  -0.20
      2 C2  1.5 0.0 0.0 C.2  1 LIG   0.62
      3 O1  2.1 1.0 0.0 O.2  1 LIG  -0.55
      4 N1  2.1 -1.2 0.0 N.am 1 LIG  -0.70
@<TRIPOS>BOND
     1 1 2 1
     2 2 3 2
     3 2 4 am
"""
    run("MOL2 'am' bond is single, not aromatic", lambda: str(sdfrust.parse_mol2_string(mol2_am).bonds[2].order), str(sdfrust.parse_mol2_string(mol2_am).bonds[0].order))
    run("MOL2 partial charges not rounded into formal charges", lambda: [a.formal_charge for a in sdfrust.parse_mol2_string(mol2_am).atoms], [0, 0, 0, 0])

    def chg_roundtrip():
        m = sdfrust.parse_sdf_string(v2000)
        atoms = m.atoms
        m2 = sdfrust.Molecule(m.name)
        for i, a in enumerate(atoms):
            na = sdfrust.Atom(i, a.element, a.x, a.y, a.z)
            na.formal_charge = {0: 1, 1: 4}.get(i, 0)
            m2.add_atom(na)
        block = sdfrust.write_sdf_string(m2)
        rd = Chem.MolFromMolBlock(block, sanitize=False)
        return [a.GetFormalCharge() for a in rd.GetAtoms()][:2]

    run("Charges +1 and +4 survive write -> RDKit read", chg_roundtrip, [1, 4])

    def iso_roundtrip(mass):
        rd = Chem.MolFromSmiles(f"[{mass}C]")
        block = Chem.MolToMolBlock(rd)
        ours = sdfrust.parse_sdf_string(block)
        back = Chem.MolFromMolBlock(sdfrust.write_sdf_string(ours), sanitize=False)
        return back.GetAtomWithIdx(0).GetIsotope()

    run("M ISO 13 round-trips (RDKit -> sdfrust -> RDKit)", lambda: iso_roundtrip(13), 13)
    run("M ISO 235 round-trips", lambda: iso_roundtrip(235), 235)

    def rad_roundtrip():
        block = Chem.MolToMolBlock(Chem.MolFromSmiles("[CH3]"))
        back = Chem.MolFromMolBlock(sdfrust.write_sdf_string(sdfrust.parse_sdf_string(block)), sanitize=False)
        return back.GetAtomWithIdx(0).GetNumRadicalElectrons()

    run("M RAD doublet round-trips", rad_roundtrip, 1)

    def big_v2000():
        m = sdfrust.Molecule("big")
        for i in range(1000):
            m.add_atom(sdfrust.Atom(i, "C", float(i), 0.0, 0.0))
        try:
            s = sdfrust.write_sdf_string(m)
        except BaseException as exc:
            return "error" if not isinstance(exc, KeyboardInterrupt) else None
        return "valid" if Chem.MolFromMolBlock(s, sanitize=False) is not None and Chem.MolFromMolBlock(s, sanitize=False).GetNumAtoms() == 1000 else "corrupt"

    run("V2000 writer with 1000 atoms errors or writes valid file", big_v2000, "error")

    run("ecfp(n_bits=0) raises a clean error", lambda: _raises_clean(lambda: sdfrust.parse_sdf_string(v2000).ecfp(2, 0)), True)
    run("add_bond with out-of-range index raises", lambda: _raises_clean(lambda: _bad_bond()), True)

    if sdfrust.gzip_enabled():
        def multi_member():
            p = workdir / "multi.sdf.gz"
            p.write_bytes(gzip.compress(v2000.encode()) + gzip.compress(v2000.encode()))
            return len(sdfrust.parse_sdf_file_multi(str(p)))

        run("Multi-member gzip reads all members", multi_member, 2)

    def prop_order():
        m = sdfrust.parse_sdf_string(v2000)
        for k in ["zeta", "alpha", "mid", "beta"]:
            m.set_property(k, k)
        outs = {sdfrust.write_sdf_string(m) for _ in range(5)}
        return len(outs)

    run("Property output is deterministic", prop_order, 1)
    return results


def _bad_bond():
    m = sdfrust.Molecule("x")
    m.add_atom(sdfrust.Atom(0, "C", 0.0, 0.0, 0.0))
    m.add_bond(sdfrust.Bond(0, 7, sdfrust.BondOrder.single()))
    sdfrust.write_sdf_v3000_string(m)
    sdfrust.write_sdf_string(m)


def _raises_clean(fn):
    try:
        fn()
    except Exception:
        return True
    except BaseException as exc:  # pyo3 PanicException
        if isinstance(exc, KeyboardInterrupt):
            raise
        return f"PANIC: {exc}"[:120]
    return False


def _iter_count(workdir, text):
    p = workdir / "iter.sdf"
    p.write_text(text)
    n = 0
    it = sdfrust.iter_sdf_file(str(p))
    while True:
        try:
            next(it)
            n += 1
        except StopIteration:
            return n
        except Exception:
            continue


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "sdfrust-python/benchmarks/reports"))
    ap.add_argument("--label", default="comparison")
    ap.add_argument("--pdbbind", type=int, default=0, help="also compare N PDBbind ligands")
    ap.add_argument("--smiles", action="append", default=[],
                    help="SMILES file to convert into an extra corpus (repeatable)")
    ap.add_argument("--smiles-limit", type=int, default=2000)
    ap.add_argument("--extra", action="append", default=[], help="extra SDF file (repeatable)")
    ap.add_argument("--rdkit-data", action="store_true",
                    help="add RDKit's bundled NCI first_5K.smi and EGFR/CDK2 ligand SDFs")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t = Tally()
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        smiles, extra = list(args.smiles), list(args.extra)
        if args.rdkit_data:
            from rdkit import RDConfig
            smiles.append(os.path.join(RDConfig.RDDataDir, "NCI", "first_5K.smi"))
            extra += [os.path.join(RDConfig.RDContribDir, "PBF", "testData", "egfr.sdf"),
                      os.path.join(RDConfig.RDContribDir, "Fastcluster", "testdata", "cdk2.sdf")]
        files = build_corpus(work, args.pdbbind, smiles, args.smiles_limit, extra)
        fps = []
        for path in files:
            try:
                ours_list = sdfrust.parse_auto_file_multi(str(path))
            except BaseException as exc:
                if isinstance(exc, KeyboardInterrupt):
                    raise
                t.crash("parse.file", path.name, exc)
                continue
            rd_list = list(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True))
            t.check("parse.record_count", len(ours_list), len(rd_list), path.name)
            ob = obabel_props(path) if len(rd_list) == 1 else None
            raw_list = list(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False))
            for i, (ours, rd) in enumerate(zip(ours_list, rd_list)):
                compare_molecule(t, path, i, ours, rd, ob, raw_list[i] if i < len(raw_list) else None)
                if rd is not None and "__noh" not in path.name:
                    try:
                        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                        fps.append((ours.ecfp(2, 2048), gen.GetFingerprint(rd)))
                    except BaseException as exc:
                        if isinstance(exc, KeyboardInterrupt):
                            raise
        edges = edge_checks(work)

    # ECFP: does sdfrust similarity rank pairs the way RDKit Morgan does?
    ecfp = {}
    if len(fps) > 3:
        if len(fps) > 400:  # pairwise cost is quadratic; a fixed sample is enough
            import random
            fps = random.Random(0).sample(fps, 400)
        ours_s, ref_s = [], []
        for i in range(len(fps)):
            for j in range(i + 1, len(fps)):
                ours_s.append(tanimoto(fps[i][0], fps[j][0]))
                ref_s.append(DataStructs.TanimotoSimilarity(fps[i][1], fps[j][1]))
        ecfp = {"pairs": len(ours_s), "spearman_vs_rdkit_morgan": round(spearman(ours_s, ref_s), 3)}

    features = sorted(t.total)
    summary = {
        f: {
            "agree": t.ok[f],
            "total": t.total[f],
            "pct": round(100 * t.ok[f] / t.total[f], 1) if t.total[f] else None,
            "max_abs_err": round(max(t.numeric[f]), 4) if t.numeric[f] else None,
            "mae": round(sum(t.numeric[f]) / len(t.numeric[f]), 4) if t.numeric[f] else None,
        }
        for f in features
    }
    report = {"label": args.label, "files": len(files), "summary": summary, "ecfp": ecfp,
              "edge_checks": edges, "worst": t.worst, "diffs": {f: t.diffs[f] for f in features if t.diffs[f]},
              "references": {"rdkit": Chem.rdBase.rdkitVersion, "openbabel": bool(OBABEL),
                             "ogb": atom_to_feature_vector is not None}}
    (out / f"{args.label}.json").write_text(json.dumps(report, indent=1, default=str))

    lines = [f"# sdfrust comparison: {args.label}", "",
             f"Files: {len(files)} · RDKit {Chem.rdBase.rdkitVersion} · OpenBabel: {bool(OBABEL)} · ogb: {atom_to_feature_vector is not None}", "",
             "## Edge-case checks", "", "| Check | Result | Got | Expected |", "|---|---|---|---|"]
    for e in edges:
        lines.append(f"| {e['check']} | {'PASS' if e['ok'] else 'FAIL'} | `{e['got']}` | `{e['expected']}` |")
    lines += ["", "## Feature agreement with reference toolkits", "",
              "| Feature | Agree | Total | % | MAE | Max err |", "|---|---|---|---|---|---|"]
    for f in features:
        s = summary[f]
        lines.append(f"| {f} | {s['agree']} | {s['total']} | {s['pct']} | {s['mae'] or ''} | {s['max_abs_err'] or ''} |")
    if ecfp:
        lines += ["", f"ECFP pairwise-similarity Spearman vs RDKit Morgan (r=2, 2048): **{ecfp['spearman_vs_rdkit_morgan']}** over {ecfp['pairs']} pairs"]
    if t.worst:
        lines += ["", "## Largest numeric errors", ""]
        for f, rows in sorted(t.worst.items()):
            lines.append(f"- **{f}**: " + "; ".join(f"`{c}` {o} vs {r}" for _, c, o, r in rows[:3]))
    lines += ["", "## Sample disagreements", ""]
    for f in features:
        if t.diffs[f]:
            lines.append(f"### {f}")
            for d in t.diffs[f][:8]:
                lines.append(f"- `{d['ctx']}`: sdfrust=`{d['sdfrust']}` reference=`{d['reference']}`")
            lines.append("")
    (out / f"{args.label}.md").write_text("\n".join(lines))
    print("\n".join(lines[: 12 + len(edges) + len(features) + 6]))


if __name__ == "__main__":
    sys.exit(main())
