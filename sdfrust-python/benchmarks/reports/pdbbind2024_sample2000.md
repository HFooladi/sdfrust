# sdfrust comparison: pdbbind_sample

Files: 2136 · RDKit 2025.09.5 · OpenBabel: True · ogb: True

## Edge-case checks

| Check | Result | Got | Expected |
|---|---|---|---|
| V3000 via parse_sdf_string gives atoms (or error), not empty | PASS | `6` | `6` |
| Mixed V2000+V3000 file parses both records | PASS | `[21, 6]` | `[21, 6]` |
| Trailing blank line after $$$$ | PASS | `1` | `1` |
| Iterator recovers after malformed record | PASS | `2` | `2` |
| Non-ASCII junk raises an error, not a panic | PASS | `True` | `True` |
| Huge V3000 count raises an error, not an abort | PASS | `True` | `True` |
| MOL2 'am' bond is single, not aromatic | PASS | `'single'` | `'single'` |
| MOL2 partial charges not rounded into formal charges | PASS | `[0, 0, 0, 0]` | `[0, 0, 0, 0]` |
| Charges +1 and +4 survive write -> RDKit read | PASS | `[1, 4]` | `[1, 4]` |
| M ISO 13 round-trips (RDKit -> sdfrust -> RDKit) | PASS | `13` | `13` |
| M ISO 235 round-trips | PASS | `235` | `235` |
| M RAD doublet round-trips | PASS | `1` | `1` |
| V2000 writer with 1000 atoms errors or writes valid file | PASS | `'error'` | `'error'` |
| ecfp(n_bits=0) raises a clean error | PASS | `True` | `True` |
| add_bond with out-of-range index raises | PASS | `True` | `True` |
| Multi-member gzip reads all members | PASS | `2` | `2` |
| Property output is deterministic | PASS | `1` | `1` |

## Feature agreement with reference toolkits

| Feature | Agree | Total | % | MAE | Max err |
|---|---|---|---|---|---|
| atom.aromatic | 122731 | 122731 | 100.0 |  |  |
| atom.chiral_tag | 122631 | 122731 | 99.9 |  |  |
| atom.hybridization | 122728 | 122731 | 100.0 |  |  |
| atom.implicit_h | 122731 | 122731 | 100.0 |  |  |
| atom.in_ring | 122731 | 122731 | 100.0 |  |  |
| atom.total_degree | 122731 | 122731 | 100.0 |  |  |
| bond.aromatic | 125649 | 125649 | 100.0 |  |  |
| bond.conjugated | 125648 | 125649 | 100.0 |  |  |
| bond.in_ring | 125649 | 125649 | 100.0 |  |  |
| charge.gasteiger_vs_openbabel | 96968 | 120301 | 80.6 | 0.0379 | 0.7854 |
| charge.gasteiger_vs_rdkit | 122582 | 122611 | 100.0 |  | 0.4522 |
| desc.exact_mass | 2140 | 2140 | 100.0 |  |  |
| desc.formula | 2140 | 2140 | 100.0 |  |  |
| desc.molecular_weight | 2140 | 2140 | 100.0 |  |  |
| desc.ring_count | 2140 | 2140 | 100.0 |  |  |
| desc.rotatable_bonds | 2140 | 2140 | 100.0 |  |  |
| desc.sssr_sizes | 2140 | 2140 | 100.0 |  |  |
| ecfp.morgan_bits_r2 | 2138 | 2140 | 99.9 |  |  |
| ecfp.morgan_bits_r3 | 2138 | 2140 | 99.9 |  |  |
| ogb.atom.aromatic | 122731 | 122731 | 100.0 |  |  |
| ogb.atom.atomic_num | 122731 | 122731 | 100.0 |  |  |
| ogb.atom.chirality | 122631 | 122731 | 99.9 |  |  |
| ogb.atom.degree | 122731 | 122731 | 100.0 |  |  |
| ogb.atom.formal_charge | 122727 | 122731 | 100.0 |  |  |
| ogb.atom.hybridization | 122728 | 122731 | 100.0 |  |  |
| ogb.atom.in_ring | 122731 | 122731 | 100.0 |  |  |
| ogb.atom.num_h | 122731 | 122731 | 100.0 |  |  |
| ogb.atom.radical_e | 122731 | 122731 | 100.0 |  |  |
| ogb.bond.bond_type | 125647 | 125649 | 100.0 |  |  |
| ogb.bond.conjugated | 125648 | 125649 | 100.0 |  |  |
| ogb.bond.stereo | 125648 | 125649 | 100.0 |  |  |
| openbabel.formula | 1587 | 2134 | 74.4 |  |  |
| openbabel.molecular_weight | 2134 | 2134 | 100.0 | 0.007 | 0.045 |
| parse.atom_count | 2140 | 2140 | 100.0 |  |  |
| parse.bond_count | 2140 | 2140 | 100.0 |  |  |
| parse.formal_charges | 2138 | 2140 | 99.9 |  |  |
| parse.isotopes | 2140 | 2140 | 100.0 |  |  |
| parse.radical_electrons | 2140 | 2140 | 100.0 |  |  |
| parse.record_count | 2136 | 2136 | 100.0 |  |  |
| rdkit_can_parse | 0 | 1 | 0.0 |  |  |
| roundtrip.atom_count | 2140 | 2140 | 100.0 |  |  |
| roundtrip.bond_orders | 2140 | 2140 | 100.0 |  |  |
| roundtrip.formal_charges | 2138 | 2140 | 99.9 |  |  |
| roundtrip.isotopes | 2140 | 2140 | 100.0 |  |  |
| roundtrip.radicals | 2140 | 2140 | 100.0 |  |  |

ECFP pairwise-similarity Spearman vs RDKit Morgan (r=2, 2048): **1.0** over 79800 pairs

## Largest numeric errors

- **charge.gasteiger_vs_openbabel**: `1sqc_ligand.sdf#0 atom 1` -0.3669 vs 0.4185; `v3000_charged.sdf#0 atom 0` -0.3689 vs 0.377; `1gvx_ligand.sdf#0 atom 0` -0.336 vs 0.4078
- **charge.gasteiger_vs_rdkit**: `1mue_ligand.sdf#0 atom 4` -0.1662 vs -0.6183; `1ksn_ligand.sdf#0 atom 0` -0.1667 vs -0.6187; `1epq_ligand.sdf#0 atom 1` -0.1509 vs 0.0141
- **desc.exact_mass**: `1bq3_ligand.sdf#0` 648.7813 vs 648.7813; `1e03_ligand.sdf#0` 1648.8584 vs 1648.8584; `1fao_ligand.sdf#0` 491.8705 vs 491.8705
- **desc.molecular_weight**: `z_2_butene__noh.sdf#0` 56.108 vs 56.108; `z_2_butene__h.sdf#0` 56.108 vs 56.108; `v3000_stereo.sdf#0` 147.374 vs 147.374
- **openbabel.molecular_weight**: `1kat_ligand.sdf#0` 2339.605 vs 2339.56; `1dva_ligand.sdf#0` 2154.421 vs 2154.38; `1lxh_ligand.sdf#0` 2339.67 vs 2339.63

## Sample disagreements

### atom.chiral_tag
- `adamantane__h.sdf#0 atom 1 C`: sdfrust=`Unspecified` reference=`CCW`
- `adamantane__h.sdf#0 atom 3 C`: sdfrust=`Unspecified` reference=`CCW`
- `adamantane__h.sdf#0 atom 5 C`: sdfrust=`Unspecified` reference=`CW`
- `adamantane__h.sdf#0 atom 7 C`: sdfrust=`Unspecified` reference=`CCW`
- `cubane__h.sdf#0 atom 0 C`: sdfrust=`Unspecified` reference=`CCW`
- `cubane__h.sdf#0 atom 1 C`: sdfrust=`Unspecified` reference=`CW`
- `cubane__h.sdf#0 atom 2 C`: sdfrust=`Unspecified` reference=`CW`
- `cubane__h.sdf#0 atom 3 C`: sdfrust=`Unspecified` reference=`CCW`

### atom.hybridization
- `1epq_ligand.sdf#0 atom 1 S`: sdfrust=`SP3D` reference=`SP3`
- `1yyy_ligand.sdf#0 atom 0 S`: sdfrust=`SP3D` reference=`SP3`
- `1zzz_ligand.sdf#0 atom 0 S`: sdfrust=`SP3D` reference=`SP3`

### bond.conjugated
- `2bow_ligand.sdf#0 bond 3 0-24`: sdfrust=`True` reference=`False`

### charge.gasteiger_vs_openbabel
- `acetaminophen.sdf#0 atom 2`: sdfrust=`-0.3263` reference=`-0.2857`
- `v3000_benzene.sdf#0 atom 0`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_benzene.sdf#0 atom 1`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_benzene.sdf#0 atom 2`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_benzene.sdf#0 atom 3`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_benzene.sdf#0 atom 4`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_benzene.sdf#0 atom 5`: sdfrust=`-0.0623` reference=`0.0`
- `v3000_charged.sdf#0 atom 0`: sdfrust=`-0.3689` reference=`0.377`

### charge.gasteiger_vs_rdkit
- `1epq_ligand.sdf#0 atom 1`: sdfrust=`-0.1509` reference=`0.0141`
- `1epq_ligand.sdf#0 atom 2`: sdfrust=`-0.288` reference=`-0.3189`
- `1epq_ligand.sdf#0 atom 3`: sdfrust=`-0.288` reference=`-0.3189`
- `1epq_ligand.sdf#0 atom 8`: sdfrust=`-0.2095` reference=`-0.23`
- `1yyy_ligand.sdf#0 atom 0`: sdfrust=`-0.1509` reference=`0.0141`
- `1yyy_ligand.sdf#0 atom 1`: sdfrust=`-0.288` reference=`-0.3189`
- `1yyy_ligand.sdf#0 atom 2`: sdfrust=`-0.288` reference=`-0.3189`
- `1yyy_ligand.sdf#0 atom 8`: sdfrust=`-0.2098` reference=`-0.2303`

### ecfp.morgan_bits_r2
- `1ksn_ligand.sdf#0`: sdfrust=`[1, 24, 57, 78, 80, 91, 110, 147, 174, 231, 315, 355, 363, 372, 401, 430, 436, 457, 480, 486, 537, 544, 559, 650, 654, 695, 807, 935, 946, 962, 1037, 1057, 1073, 1083, 1148, 1152, 1160, 1171, 1222, 1230, 1231, 1250, 1261, 1349, 1359, 1380, 1385, 1398, 1488, 1562, 1564, 1652, 1673, 1680, 1700, 1754, 1757, 1763, 1789, 1866, 1869, 1873, 1917, 1932, 1961]` reference=`[1, 24, 78, 80, 91, 110, 147, 174, 231, 315, 355, 363, 372, 398, 430, 436, 457, 463, 486, 537, 559, 650, 654, 695, 715, 791, 807, 946, 962, 1037, 1057, 1073, 1083, 1148, 1152, 1160, 1171, 1222, 1230, 1231, 1250, 1261, 1314, 1349, 1359, 1380, 1385, 1398, 1488, 1514, 1562, 1652, 1673, 1680, 1700, 1754, 1757, 1763, 1786, 1795, 1866, 1869, 1873, 1917, 1932, 1961]`
- `1mue_ligand.sdf#0`: sdfrust=`[5, 41, 80, 106, 114, 212, 314, 363, 378, 504, 523, 561, 564, 650, 699, 708, 729, 794, 802, 807, 847, 872, 935, 1060, 1073, 1075, 1099, 1148, 1152, 1165, 1166, 1168, 1276, 1321, 1349, 1380, 1381, 1445, 1452, 1453, 1463, 1558, 1564, 1582, 1599, 1652, 1656, 1683, 1755, 1789, 1806, 1808, 1840, 1855, 1873, 1885, 1892, 1917, 1928, 1932, 2011]` reference=`[5, 41, 80, 114, 212, 314, 363, 378, 398, 463, 504, 523, 540, 554, 561, 564, 650, 699, 708, 715, 794, 802, 807, 872, 918, 923, 935, 1060, 1073, 1075, 1099, 1145, 1148, 1152, 1166, 1276, 1314, 1321, 1349, 1380, 1381, 1445, 1452, 1453, 1463, 1558, 1579, 1599, 1652, 1656, 1683, 1698, 1755, 1806, 1840, 1855, 1873, 1892, 1917, 1928, 1932, 2011]`

### ecfp.morgan_bits_r3
- `1ksn_ligand.sdf#0`: sdfrust=`[1, 24, 30, 57, 78, 80, 91, 110, 147, 174, 231, 315, 325, 355, 363, 372, 377, 396, 401, 406, 430, 436, 457, 480, 486, 537, 544, 556, 559, 567, 582, 650, 654, 695, 697, 807, 919, 935, 946, 962, 1037, 1057, 1073, 1083, 1140, 1148, 1152, 1160, 1171, 1218, 1222, 1230, 1231, 1250, 1261, 1349, 1359, 1380, 1385, 1398, 1449, 1477, 1488, 1562, 1564, 1567, 1578, 1652, 1673, 1680, 1700, 1748, 1754, 1757, 1763, 1789, 1866, 1869, 1873, 1887, 1892, 1917, 1932, 1947, 1961]` reference=`[1, 24, 30, 78, 80, 91, 110, 147, 174, 231, 315, 325, 355, 363, 372, 377, 396, 398, 430, 436, 457, 463, 486, 537, 556, 559, 567, 641, 650, 654, 695, 697, 715, 791, 807, 919, 946, 962, 1037, 1057, 1073, 1078, 1083, 1140, 1148, 1152, 1160, 1171, 1218, 1222, 1230, 1231, 1250, 1261, 1314, 1349, 1359, 1380, 1385, 1398, 1449, 1488, 1514, 1562, 1567, 1578, 1595, 1652, 1673, 1680, 1700, 1748, 1754, 1757, 1763, 1786, 1795, 1838, 1866, 1869, 1873, 1887, 1892, 1917, 1932, 1961]`
- `1mue_ligand.sdf#0`: sdfrust=`[5, 41, 80, 95, 106, 114, 134, 212, 314, 363, 378, 411, 476, 486, 504, 523, 561, 564, 650, 682, 692, 699, 708, 729, 794, 802, 807, 847, 872, 935, 961, 965, 1060, 1073, 1075, 1099, 1113, 1116, 1148, 1152, 1165, 1166, 1168, 1181, 1191, 1192, 1276, 1310, 1315, 1321, 1349, 1380, 1381, 1445, 1452, 1453, 1463, 1558, 1564, 1565, 1582, 1591, 1599, 1601, 1652, 1656, 1683, 1724, 1755, 1784, 1789, 1806, 1808, 1840, 1855, 1873, 1885, 1892, 1917, 1928, 1932, 1980, 2011]` reference=`[5, 9, 41, 80, 95, 114, 212, 314, 363, 373, 378, 389, 398, 411, 463, 476, 486, 504, 523, 540, 554, 561, 564, 650, 682, 692, 699, 708, 715, 794, 802, 807, 872, 915, 918, 923, 935, 965, 1060, 1073, 1075, 1099, 1113, 1145, 1148, 1152, 1166, 1191, 1276, 1310, 1314, 1315, 1321, 1349, 1380, 1381, 1445, 1452, 1453, 1463, 1486, 1558, 1579, 1599, 1601, 1652, 1656, 1683, 1698, 1713, 1724, 1741, 1755, 1806, 1812, 1840, 1855, 1873, 1892, 1917, 1928, 1932, 1980, 2011]`

### ogb.atom.chirality
- `adamantane__h.sdf#0 atom 1 C`: sdfrust=`0` reference=`2`
- `adamantane__h.sdf#0 atom 3 C`: sdfrust=`0` reference=`2`
- `adamantane__h.sdf#0 atom 5 C`: sdfrust=`0` reference=`1`
- `adamantane__h.sdf#0 atom 7 C`: sdfrust=`0` reference=`2`
- `cubane__h.sdf#0 atom 0 C`: sdfrust=`0` reference=`2`
- `cubane__h.sdf#0 atom 1 C`: sdfrust=`0` reference=`1`
- `cubane__h.sdf#0 atom 2 C`: sdfrust=`0` reference=`1`
- `cubane__h.sdf#0 atom 3 C`: sdfrust=`0` reference=`2`

### ogb.atom.formal_charge
- `1ksn_ligand.sdf#0 atom 0 O`: sdfrust=`5` reference=`4`
- `1ksn_ligand.sdf#0 atom 7 N`: sdfrust=`5` reference=`6`
- `1mue_ligand.sdf#0 atom 4 O`: sdfrust=`5` reference=`4`
- `1mue_ligand.sdf#0 atom 10 N`: sdfrust=`5` reference=`6`

### ogb.atom.hybridization
- `1epq_ligand.sdf#0 atom 1 S`: sdfrust=`3` reference=`2`
- `1yyy_ligand.sdf#0 atom 0 S`: sdfrust=`3` reference=`2`
- `1zzz_ligand.sdf#0 atom 0 S`: sdfrust=`3` reference=`2`

### ogb.bond.bond_type
- `1ksn_ligand.sdf#0 bond 0`: sdfrust=`1` reference=`0`
- `1mue_ligand.sdf#0 bond 4`: sdfrust=`1` reference=`0`

### ogb.bond.conjugated
- `2bow_ligand.sdf#0 bond 3`: sdfrust=`1` reference=`0`

### ogb.bond.stereo
- `1fd0_ligand.sdf#0 bond 4`: sdfrust=`2` reference=`1`

### openbabel.formula
- `v3000_benzene.sdf#0`: sdfrust=`C6H6` reference=`C6`
- `cisplatin__h.sdf#0`: sdfrust=`H6Cl2N2Pt` reference=`Cl2H6N2Pt`
- `cisplatin__noh.sdf#0`: sdfrust=`H6Cl2N2Pt` reference=`Cl2H6N2Pt`
- `deuterated_benzene__h.sdf#0`: sdfrust=`C6H6` reference=`C6H5D`
- `deuterated_benzene__noh.sdf#0`: sdfrust=`C6H6` reference=`C6H5D`
- `ferrocene_like__h.sdf#0`: sdfrust=`Fe+2` reference=`Fe++`
- `ferrocene_like__noh.sdf#0`: sdfrust=`Fe+2` reference=`Fe++`
- `11gs_ligand.sdf#0`: sdfrust=`C23H27Cl2N3O10S-2` reference=`C23H27Cl2N3O10S--`

### parse.formal_charges
- `1ksn_ligand.sdf#0`: sdfrust=`[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]` reference=`[-1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`
- `1mue_ligand.sdf#0`: sdfrust=`[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]` reference=`[0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`

### rdkit_can_parse
- `1lvk_ligand.sdf#0`: sdfrust=`True` reference=`False`

### roundtrip.formal_charges
- `1ksn_ligand.sdf#0`: sdfrust=`[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]` reference=`[-1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`
- `1mue_ligand.sdf#0`: sdfrust=`[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]` reference=`[0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]`
