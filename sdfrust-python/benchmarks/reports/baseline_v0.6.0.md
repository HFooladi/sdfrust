# sdfrust comparison: baseline_v0.6.0

Files: 136 · RDKit 2025.09.5 · OpenBabel: True · ogb: True

## Edge-case checks

| Check | Result | Got | Expected |
|---|---|---|---|
| V3000 via parse_sdf_string gives atoms (or error), not empty | FAIL | `0` | `12` |
| Mixed V2000+V3000 file parses both records | FAIL | `[21, 0]` | `[21, 12]` |
| Trailing blank line after $$$$ | FAIL | `ValueError: Missing required section: header` | `1` |
| Iterator recovers after malformed record | FAIL | `0` | `2` |
| Non-ASCII junk raises an error, not a panic | FAIL | `"PANIC: end byte index 3 is not a char boundary; it is inside 'Å' (bytes 2..4 of string)"` | `True` |
| Huge V3000 count raises an error, not an abort | PASS | `True` | `True` |
| MOL2 'am' bond is single, not aromatic | FAIL | `'aromatic'` | `'single'` |
| MOL2 partial charges not rounded into formal charges | FAIL | `[0, 1, -1, -1]` | `[0, 0, 0, 0]` |
| Charges +1 and +4 survive write -> RDKit read | FAIL | `[0, 4]` | `[1, 4]` |
| M ISO 13 round-trips (RDKit -> sdfrust -> RDKit) | PASS | `13` | `13` |
| M ISO 235 round-trips | FAIL | `0` | `235` |
| M RAD doublet round-trips | FAIL | `0` | `1` |
| V2000 writer with 1000 atoms errors or writes valid file | FAIL | `'corrupt'` | `'error'` |
| ecfp(n_bits=0) raises a clean error | FAIL | `'PANIC: attempt to calculate the remainder with a divisor of zero'` | `True` |
| add_bond with out-of-range index raises | FAIL | `'PANIC: index out of bounds: the len is 1 but the index is 7'` | `True` |
| Multi-member gzip reads all members | FAIL | `1` | `2` |
| Property output is deterministic | PASS | `1` | `1` |

## Feature agreement with reference toolkits

| Feature | Agree | Total | % | MAE | Max err |
|---|---|---|---|---|---|
| atom.aromatic | 1995 | 2059 | 96.9 |  |  |
| atom.chiral_tag | 2043 | 2059 | 99.2 |  |  |
| atom.hybridization | 2021 | 2059 | 98.2 |  |  |
| atom.implicit_h | 2057 | 2059 | 99.9 |  |  |
| atom.in_ring | 2059 | 2059 | 100.0 |  |  |
| atom.total_degree | 2057 | 2059 | 99.9 |  |  |
| bond.aromatic | 1995 | 2075 | 96.1 |  |  |
| bond.conjugated | 2041 | 2075 | 98.4 |  |  |
| bond.in_ring | 2075 | 2075 | 100.0 |  |  |
| charge.gasteiger_vs_openbabel | 1029 | 1794 | 57.4 | 0.0333 | 1.1107 |
| charge.gasteiger_vs_rdkit | 897 | 2015 | 44.5 | 0.0388 | 2.0 |
| desc.exact_mass | 77 | 141 | 54.6 | 8.589 | 349.0343 |
| desc.formula | 73 | 141 | 51.8 |  |  |
| desc.molecular_weight | 77 | 141 | 54.6 | 8.5899 | 349.0319 |
| desc.ring_count | 141 | 141 | 100.0 |  |  |
| desc.rotatable_bonds | 121 | 141 | 85.8 |  |  |
| desc.sssr_sizes | 141 | 141 | 100.0 |  |  |
| ogb.atom.aromatic | 1995 | 2059 | 96.9 |  |  |
| ogb.atom.atomic_num | 0 | 2059 | 0.0 |  |  |
| ogb.atom.chirality | 2043 | 2059 | 99.2 |  |  |
| ogb.atom.degree | 1706 | 2059 | 82.9 |  |  |
| ogb.atom.formal_charge | 2059 | 2059 | 100.0 |  |  |
| ogb.atom.hybridization | 18 | 2059 | 0.9 |  |  |
| ogb.atom.in_ring | 2059 | 2059 | 100.0 |  |  |
| ogb.atom.num_h | 2057 | 2059 | 99.9 |  |  |
| ogb.atom.radical_e | 2057 | 2059 | 99.9 |  |  |
| ogb.bond.bond_type | 1995 | 2075 | 96.1 |  |  |
| ogb.bond.conjugated | 2041 | 2075 | 98.4 |  |  |
| ogb.bond.stereo | 2071 | 2075 | 99.8 |  |  |
| openbabel.formula | 66 | 135 | 48.9 |  |  |
| openbabel.molecular_weight | 71 | 135 | 52.6 | 8.9716 | 349.034 |
| parse.atom_count | 141 | 141 | 100.0 |  |  |
| parse.bond_count | 141 | 141 | 100.0 |  |  |
| parse.formal_charges | 141 | 141 | 100.0 |  |  |
| parse.isotopes | 139 | 141 | 98.6 |  |  |
| parse.record_count | 136 | 136 | 100.0 |  |  |
| roundtrip.atom_count | 141 | 141 | 100.0 |  |  |
| roundtrip.bond_orders | 141 | 141 | 100.0 |  |  |
| roundtrip.formal_charges | 141 | 141 | 100.0 |  |  |
| roundtrip.isotopes | 139 | 141 | 98.6 |  |  |
| roundtrip.radicals | 139 | 141 | 98.6 |  |  |

ECFP pairwise-similarity Spearman vs RDKit Morgan (r=2, 2048): **0.915** over 3160 pairs

## Largest numeric errors

- **charge.gasteiger_vs_openbabel**: `v3000_charged.sdf#0 atom 0` -0.7337 vs 0.377; `glycine_zw__h.sdf#0 atom 0` -0.6263 vs 0.3852; `tetramethylammonium__h.sdf#0 atom 1` -0.3713 vs 0.3725
- **charge.gasteiger_vs_rdkit**: `ferrocene_like__noh.sdf#0 atom 0` 0.0 vs 2.0; `ferrocene_like__h.sdf#0 atom 0` 0.0 vs 2.0; `glycine_zw__noh.sdf#0 atom 0` 0.4546 vs -0.3528
- **desc.exact_mass**: `uranium235_f6__noh.sdf#0` 0.0 vs 349.0343; `uranium235_f6__h.sdf#0` 0.0 vs 349.0343; `paclitaxel_core__noh.sdf#0` 547.9441 vs 586.2414
- **desc.molecular_weight**: `uranium235_f6__noh.sdf#0` 0.0 vs 349.0319; `uranium235_f6__h.sdf#0` 0.0 vs 349.0319; `paclitaxel_core__noh.sdf#0` 548.33 vs 586.634
- **openbabel.molecular_weight**: `uranium235_f6__noh.sdf#0` 0.0 vs 349.034; `uranium235_f6__h.sdf#0` 0.0 vs 349.034; `paclitaxel_core__noh.sdf#0` 548.33 vs 586.627

## Sample disagreements

### atom.aromatic
- `anthracene__h.sdf#0 atom 4 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 5 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 6 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 7 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 8 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 9 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 10 C`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 atom 11 C`: sdfrust=`False` reference=`True`

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
- `allene__h.sdf#0 atom 1 C`: sdfrust=`SP2` reference=`SP`
- `allene__noh.sdf#0 atom 1 C`: sdfrust=`SP2` reference=`SP`
- `anilinium__h.sdf#0 atom 0 N`: sdfrust=`SP2` reference=`SP3`
- `anilinium__noh.sdf#0 atom 0 N`: sdfrust=`SP2` reference=`SP3`
- `azide__h.sdf#0 atom 2 N`: sdfrust=`SP2` reference=`SP`
- `azide__noh.sdf#0 atom 2 N`: sdfrust=`SP2` reference=`SP`
- `boronic_acid__h.sdf#0 atom 1 B`: sdfrust=`SP3` reference=`SP2`
- `boronic_acid__noh.sdf#0 atom 1 B`: sdfrust=`SP3` reference=`SP2`

### atom.implicit_h
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`1` reference=`0`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`

### atom.total_degree
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`

### bond.aromatic
- `anthracene__h.sdf#0 bond 3 3-4`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 4 4-5`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 5 5-6`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 6 6-7`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 7 7-8`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 8 8-9`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 9 9-10`: sdfrust=`False` reference=`True`
- `anthracene__h.sdf#0 bond 10 10-11`: sdfrust=`False` reference=`True`

### bond.conjugated
- `anilinium__h.sdf#0 bond 0 0-1`: sdfrust=`True` reference=`False`
- `anilinium__noh.sdf#0 bond 0 0-1`: sdfrust=`True` reference=`False`
- `benzonitrile__h.sdf#0 bond 1 1-2`: sdfrust=`False` reference=`True`
- `benzonitrile__noh.sdf#0 bond 1 1-2`: sdfrust=`False` reference=`True`
- `cyclopentadienide__h.sdf#0 bond 0 0-1`: sdfrust=`False` reference=`True`
- `cyclopentadienide__h.sdf#0 bond 4 4-0`: sdfrust=`False` reference=`True`
- `cyclopentadienide__noh.sdf#0 bond 0 0-1`: sdfrust=`False` reference=`True`
- `cyclopentadienide__noh.sdf#0 bond 4 4-0`: sdfrust=`False` reference=`True`

### charge.gasteiger_vs_openbabel
- `acetaminophen.sdf#0 atom 0`: sdfrust=`-0.4094` reference=`-0.5068`
- `acetaminophen.sdf#0 atom 2`: sdfrust=`-0.2566` reference=`-0.2857`
- `acetaminophen.sdf#0 atom 3`: sdfrust=`0.0068` reference=`0.0336`
- `acetaminophen.sdf#0 atom 4`: sdfrust=`-0.0605` reference=`-0.037`
- `acetaminophen.sdf#0 atom 5`: sdfrust=`-0.0605` reference=`-0.037`
- `acetaminophen.sdf#0 atom 6`: sdfrust=`-0.0508` reference=`-0.0184`
- `acetaminophen.sdf#0 atom 7`: sdfrust=`-0.0508` reference=`-0.0184`
- `acetaminophen.sdf#0 atom 8`: sdfrust=`0.0508` reference=`0.117`

### charge.gasteiger_vs_rdkit
- `acetaminophen.sdf#0 atom 0`: sdfrust=`-0.4094` reference=`-0.508`
- `acetaminophen.sdf#0 atom 2`: sdfrust=`-0.2566` reference=`-0.3263`
- `acetaminophen.sdf#0 atom 3`: sdfrust=`0.0068` reference=`0.0409`
- `acetaminophen.sdf#0 atom 4`: sdfrust=`-0.0605` reference=`-0.0349`
- `acetaminophen.sdf#0 atom 5`: sdfrust=`-0.0605` reference=`-0.0349`
- `acetaminophen.sdf#0 atom 6`: sdfrust=`-0.0508` reference=`-0.0182`
- `acetaminophen.sdf#0 atom 7`: sdfrust=`-0.0508` reference=`-0.0182`
- `acetaminophen.sdf#0 atom 8`: sdfrust=`0.0508` reference=`0.1152`

### desc.exact_mass
- `v3000_benzene.sdf#0`: sdfrust=`72.0` reference=`78.047`
- `v3000_charged.sdf#0`: sdfrust=`62.0242` reference=`63.032`
- `2-pyridone__noh.sdf#0`: sdfrust=`89.998` reference=`95.0371`
- `acetamide__noh.sdf#0`: sdfrust=`53.998` reference=`59.0371`
- `acetate__noh.sdf#0`: sdfrust=`55.9898` reference=`59.0139`
- `acetonitrile__noh.sdf#0`: sdfrust=`38.0031` reference=`41.0265`
- `adamantane__noh.sdf#0`: sdfrust=`120.0` reference=`136.1252`
- `allene__noh.sdf#0`: sdfrust=`36.0` reference=`40.0313`

### desc.formula
- `v3000_benzene.sdf#0`: sdfrust=`C6` reference=`C6H6`
- `v3000_charged.sdf#0`: sdfrust=`CH4NO2` reference=`CH5NO2`
- `2-pyridone__noh.sdf#0`: sdfrust=`C5NO` reference=`C5H5NO`
- `acetamide__noh.sdf#0`: sdfrust=`C2NO` reference=`C2H5NO`
- `acetate__h.sdf#0`: sdfrust=`C2H3O2` reference=`C2H3O2-`
- `acetate__noh.sdf#0`: sdfrust=`C2O2` reference=`C2H3O2-`
- `acetonitrile__noh.sdf#0`: sdfrust=`C2N` reference=`C2H3N`
- `adamantane__noh.sdf#0`: sdfrust=`C10` reference=`C10H16`

### desc.molecular_weight
- `v3000_benzene.sdf#0`: sdfrust=`72.066` reference=`78.114`
- `v3000_charged.sdf#0`: sdfrust=`62.048` reference=`63.056`
- `2-pyridone__noh.sdf#0`: sdfrust=`90.061` reference=`95.101`
- `acetamide__noh.sdf#0`: sdfrust=`54.028` reference=`59.068`
- `acetate__noh.sdf#0`: sdfrust=`56.02` reference=`59.044`
- `acetonitrile__noh.sdf#0`: sdfrust=`38.029` reference=`41.053`
- `adamantane__noh.sdf#0`: sdfrust=`120.11` reference=`136.238`
- `allene__noh.sdf#0`: sdfrust=`36.033` reference=`40.065`

### desc.rotatable_bonds
- `acetaminophen.sdf#0`: sdfrust=`2` reference=`1`
- `aspirin.sdf#0`: sdfrust=`3` reference=`2`
- `drug_library.sdf#3`: sdfrust=`3` reference=`2`
- `drug_library.sdf#5`: sdfrust=`2` reference=`1`
- `atorvastatin__h.sdf#0`: sdfrust=`13` reference=`12`
- `atorvastatin__noh.sdf#0`: sdfrust=`13` reference=`12`
- `benzonitrile__h.sdf#0`: sdfrust=`1` reference=`0`
- `benzonitrile__noh.sdf#0`: sdfrust=`1` reference=`0`

### ogb.atom.aromatic
- `anthracene__h.sdf#0 atom 4 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 5 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 6 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 7 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 8 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 9 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 10 C`: sdfrust=`0` reference=`1`
- `anthracene__h.sdf#0 atom 11 C`: sdfrust=`0` reference=`1`

### ogb.atom.atomic_num
- `acetaminophen.sdf#0 atom 0 O`: sdfrust=`8` reference=`7`
- `acetaminophen.sdf#0 atom 1 O`: sdfrust=`8` reference=`7`
- `acetaminophen.sdf#0 atom 2 N`: sdfrust=`7` reference=`6`
- `acetaminophen.sdf#0 atom 3 C`: sdfrust=`6` reference=`5`
- `acetaminophen.sdf#0 atom 4 C`: sdfrust=`6` reference=`5`
- `acetaminophen.sdf#0 atom 5 C`: sdfrust=`6` reference=`5`
- `acetaminophen.sdf#0 atom 6 C`: sdfrust=`6` reference=`5`
- `acetaminophen.sdf#0 atom 7 C`: sdfrust=`6` reference=`5`

### ogb.atom.chirality
- `adamantane__h.sdf#0 atom 1 C`: sdfrust=`0` reference=`2`
- `adamantane__h.sdf#0 atom 3 C`: sdfrust=`0` reference=`2`
- `adamantane__h.sdf#0 atom 5 C`: sdfrust=`0` reference=`1`
- `adamantane__h.sdf#0 atom 7 C`: sdfrust=`0` reference=`2`
- `cubane__h.sdf#0 atom 0 C`: sdfrust=`0` reference=`2`
- `cubane__h.sdf#0 atom 1 C`: sdfrust=`0` reference=`1`
- `cubane__h.sdf#0 atom 2 C`: sdfrust=`0` reference=`1`
- `cubane__h.sdf#0 atom 3 C`: sdfrust=`0` reference=`2`

### ogb.atom.degree
- `v3000_benzene.sdf#0 atom 0 C`: sdfrust=`2` reference=`3`
- `v3000_benzene.sdf#0 atom 1 C`: sdfrust=`2` reference=`3`
- `v3000_benzene.sdf#0 atom 2 C`: sdfrust=`2` reference=`3`
- `v3000_benzene.sdf#0 atom 3 C`: sdfrust=`2` reference=`3`
- `v3000_benzene.sdf#0 atom 4 C`: sdfrust=`2` reference=`3`
- `v3000_benzene.sdf#0 atom 5 C`: sdfrust=`2` reference=`3`
- `v3000_charged.sdf#0 atom 5 C`: sdfrust=`2` reference=`3`
- `2-pyridone__noh.sdf#0 atom 2 C`: sdfrust=`2` reference=`3`

### ogb.atom.hybridization
- `acetaminophen.sdf#0 atom 0 O`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 1 O`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 2 N`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 3 C`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 4 C`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 5 C`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 6 C`: sdfrust=`2` reference=`1`
- `acetaminophen.sdf#0 atom 7 C`: sdfrust=`2` reference=`1`

### ogb.atom.num_h
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`1` reference=`0`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`

### ogb.atom.radical_e
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`0` reference=`1`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`0` reference=`1`

### ogb.bond.bond_type
- `anthracene__h.sdf#0 bond 3`: sdfrust=`0` reference=`3`
- `anthracene__h.sdf#0 bond 4`: sdfrust=`1` reference=`3`
- `anthracene__h.sdf#0 bond 5`: sdfrust=`0` reference=`3`
- `anthracene__h.sdf#0 bond 6`: sdfrust=`1` reference=`3`
- `anthracene__h.sdf#0 bond 7`: sdfrust=`0` reference=`3`
- `anthracene__h.sdf#0 bond 8`: sdfrust=`1` reference=`3`
- `anthracene__h.sdf#0 bond 9`: sdfrust=`0` reference=`3`
- `anthracene__h.sdf#0 bond 10`: sdfrust=`1` reference=`3`

### ogb.bond.conjugated
- `anilinium__h.sdf#0 bond 0`: sdfrust=`1` reference=`0`
- `anilinium__noh.sdf#0 bond 0`: sdfrust=`1` reference=`0`
- `benzonitrile__h.sdf#0 bond 1`: sdfrust=`0` reference=`1`
- `benzonitrile__noh.sdf#0 bond 1`: sdfrust=`0` reference=`1`
- `cyclopentadienide__h.sdf#0 bond 0`: sdfrust=`0` reference=`1`
- `cyclopentadienide__h.sdf#0 bond 4`: sdfrust=`0` reference=`1`
- `cyclopentadienide__noh.sdf#0 bond 0`: sdfrust=`0` reference=`1`
- `cyclopentadienide__noh.sdf#0 bond 4`: sdfrust=`0` reference=`1`

### ogb.bond.stereo
- `e_2_butene__h.sdf#0 bond 1`: sdfrust=`0` reference=`2`
- `e_2_butene__noh.sdf#0 bond 1`: sdfrust=`0` reference=`2`
- `z_2_butene__h.sdf#0 bond 1`: sdfrust=`0` reference=`1`
- `z_2_butene__noh.sdf#0 bond 1`: sdfrust=`0` reference=`1`

### openbabel.formula
- `v3000_charged.sdf#0`: sdfrust=`CH4NO2` reference=`CH5NO2`
- `2-pyridone__noh.sdf#0`: sdfrust=`C5NO` reference=`C5H5NO`
- `acetamide__noh.sdf#0`: sdfrust=`C2NO` reference=`C2H5NO`
- `acetate__h.sdf#0`: sdfrust=`C2H3O2` reference=`C2H3O2-`
- `acetate__noh.sdf#0`: sdfrust=`C2O2` reference=`C2H3O2-`
- `acetonitrile__noh.sdf#0`: sdfrust=`C2N` reference=`C2H3N`
- `adamantane__noh.sdf#0`: sdfrust=`C10` reference=`C10H16`
- `allene__noh.sdf#0`: sdfrust=`C3` reference=`C3H4`

### openbabel.molecular_weight
- `v3000_benzene.sdf#0`: sdfrust=`72.066` reference=`78.11`
- `v3000_charged.sdf#0`: sdfrust=`62.048` reference=`63.0559`
- `2-pyridone__noh.sdf#0`: sdfrust=`90.061` reference=`95.0993`
- `acetamide__noh.sdf#0`: sdfrust=`54.028` reference=`59.0672`
- `acetate__noh.sdf#0`: sdfrust=`56.02` reference=`59.044`
- `acetonitrile__noh.sdf#0`: sdfrust=`38.029` reference=`41.0519`
- `adamantane__noh.sdf#0`: sdfrust=`120.11` reference=`136.234`
- `allene__noh.sdf#0`: sdfrust=`36.033` reference=`40.0639`

### parse.isotopes
- `uranium235_f6__h.sdf#0`: sdfrust=`[None, None, None, None, None, None, None]` reference=`[None, 235, None, None, None, None, None]`
- `uranium235_f6__noh.sdf#0`: sdfrust=`[None, None, None, None, None, None, None]` reference=`[None, 235, None, None, None, None, None]`

### roundtrip.isotopes
- `uranium235_f6__h.sdf#0`: sdfrust=`[0, 0, 0, 0, 0, 0, 0]` reference=`[0, 235, 0, 0, 0, 0, 0]`
- `uranium235_f6__noh.sdf#0`: sdfrust=`[0, 0, 0, 0, 0, 0, 0]` reference=`[0, 235, 0, 0, 0, 0, 0]`

### roundtrip.radicals
- `methyl_radical__h.sdf#0`: sdfrust=`[0, 0, 0, 0]` reference=`[1, 0, 0, 0]`
- `methyl_radical__noh.sdf#0`: sdfrust=`[0]` reference=`[1]`
