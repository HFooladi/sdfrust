# sdfrust comparison: phase2_baseline

Files: 140 · RDKit 2025.09.5 · OpenBabel: True · ogb: True

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
| atom.aromatic | 110779 | 113318 | 97.8 |  |  |
| atom.chiral_tag | 113089 | 113318 | 99.8 |  |  |
| atom.hybridization | 111937 | 113318 | 98.8 |  |  |
| atom.implicit_h | 113316 | 113318 | 100.0 |  |  |
| atom.in_ring | 113318 | 113318 | 100.0 |  |  |
| atom.total_degree | 113316 | 113318 | 100.0 |  |  |
| bond.aromatic | 112665 | 115958 | 97.2 |  |  |
| bond.conjugated | 114046 | 115958 | 98.4 |  |  |
| bond.in_ring | 115958 | 115958 | 100.0 |  |  |
| charge.gasteiger_vs_openbabel | 1029 | 1794 | 57.4 | 0.0333 | 1.1107 |
| charge.gasteiger_vs_rdkit | 49264 | 106524 | 46.2 | 0.0352 | 2.0 |
| desc.exact_mass | 2467 | 4553 | 54.2 | 10.9707 | 596.085 |
| desc.formula | 2340 | 4553 | 51.4 |  |  |
| desc.molecular_weight | 2467 | 4553 | 54.2 | 10.9722 | 595.927 |
| desc.ring_count | 4553 | 4553 | 100.0 |  |  |
| desc.rotatable_bonds | 3094 | 4553 | 68.0 |  |  |
| desc.sssr_sizes | 4545 | 4553 | 99.8 |  |  |
| ogb.atom.aromatic | 110779 | 113318 | 97.8 |  |  |
| ogb.atom.atomic_num | 113268 | 113318 | 100.0 |  |  |
| ogb.atom.chirality | 113089 | 113318 | 99.8 |  |  |
| ogb.atom.degree | 113316 | 113318 | 100.0 |  |  |
| ogb.atom.formal_charge | 113318 | 113318 | 100.0 |  |  |
| ogb.atom.hybridization | 111937 | 113318 | 98.8 |  |  |
| ogb.atom.in_ring | 113318 | 113318 | 100.0 |  |  |
| ogb.atom.num_h | 113316 | 113318 | 100.0 |  |  |
| ogb.atom.radical_e | 113318 | 113318 | 100.0 |  |  |
| ogb.bond.bond_type | 112665 | 115958 | 97.2 |  |  |
| ogb.bond.conjugated | 114046 | 115958 | 98.4 |  |  |
| ogb.bond.stereo | 115315 | 115958 | 99.4 |  |  |
| openbabel.formula | 66 | 135 | 48.9 |  |  |
| openbabel.molecular_weight | 71 | 135 | 52.6 | 8.9716 | 349.034 |
| parse.atom_count | 4553 | 4553 | 100.0 |  |  |
| parse.bond_count | 4553 | 4553 | 100.0 |  |  |
| parse.formal_charges | 4553 | 4553 | 100.0 |  |  |
| parse.isotopes | 4553 | 4553 | 100.0 |  |  |
| parse.radical_electrons | 4553 | 4553 | 100.0 |  |  |
| parse.record_count | 140 | 140 | 100.0 |  |  |
| roundtrip.atom_count | 4553 | 4553 | 100.0 |  |  |
| roundtrip.bond_orders | 4553 | 4553 | 100.0 |  |  |
| roundtrip.formal_charges | 4553 | 4553 | 100.0 |  |  |
| roundtrip.isotopes | 4553 | 4553 | 100.0 |  |  |
| roundtrip.radicals | 4553 | 4553 | 100.0 |  |  |

ECFP pairwise-similarity Spearman vs RDKit Morgan (r=2, 2048): **0.914** over 3103786 pairs

## Largest numeric errors

- **charge.gasteiger_vs_openbabel**: `v3000_charged.sdf#0 atom 0` -0.7337 vs 0.377; `glycine_zw__h.sdf#0 atom 0` -0.6263 vs 0.3852; `tetramethylammonium__h.sdf#0 atom 1` -0.3713 vs 0.3725
- **charge.gasteiger_vs_rdkit**: `ferrocene_like__noh.sdf#0 atom 0` 0.0 vs 2.0; `ferrocene_like__h.sdf#0 atom 0` 0.0 vs 2.0; `glycine_zw__noh.sdf#0 atom 0` 0.4546 vs -0.3528
- **desc.exact_mass**: `first_5K__h.sdf#243` 0.0 vs 596.085; `first_5K__noh.sdf#243` 0.0 vs 596.085; `first_5K__h.sdf#1243` 0.0 vs 543.0145
- **desc.molecular_weight**: `first_5K__h.sdf#243` 0.0 vs 595.927; `first_5K__noh.sdf#243` 0.0 vs 595.927; `first_5K__h.sdf#1243` 0.0 vs 545.279
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

### desc.sssr_sizes
- `first_5K__noh.sdf#862`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__noh.sdf#863`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__noh.sdf#864`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__noh.sdf#1377`: sdfrust=`[5, 5, 9]` reference=`[5, 5, 6]`
- `first_5K__h.sdf#862`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__h.sdf#863`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__h.sdf#864`: sdfrust=`[5, 6, 6, 6, 9]` reference=`[5, 6, 6, 6, 6]`
- `first_5K__h.sdf#1377`: sdfrust=`[5, 5, 9]` reference=`[5, 5, 6]`

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
- `uranium235_f6__h.sdf#0 atom 1 U`: sdfrust=`118` reference=`91`
- `uranium235_f6__noh.sdf#0 atom 1 U`: sdfrust=`118` reference=`91`
- `first_5K__noh.sdf#243 atom 4 Hg`: sdfrust=`118` reference=`79`
- `first_5K__noh.sdf#430 atom 1 Sb`: sdfrust=`118` reference=`50`
- `first_5K__noh.sdf#463 atom 6 Bi`: sdfrust=`118` reference=`82`
- `first_5K__noh.sdf#475 atom 3 Ti`: sdfrust=`118` reference=`21`
- `first_5K__noh.sdf#482 atom 1 Hg`: sdfrust=`118` reference=`79`
- `first_5K__noh.sdf#535 atom 4 Hg`: sdfrust=`118` reference=`79`

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
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`

### ogb.atom.hybridization
- `allene__h.sdf#0 atom 1 C`: sdfrust=`1` reference=`0`
- `allene__noh.sdf#0 atom 1 C`: sdfrust=`1` reference=`0`
- `anilinium__h.sdf#0 atom 0 N`: sdfrust=`1` reference=`2`
- `anilinium__noh.sdf#0 atom 0 N`: sdfrust=`1` reference=`2`
- `azide__h.sdf#0 atom 2 N`: sdfrust=`1` reference=`0`
- `azide__noh.sdf#0 atom 2 N`: sdfrust=`1` reference=`0`
- `boronic_acid__h.sdf#0 atom 1 B`: sdfrust=`2` reference=`1`
- `boronic_acid__noh.sdf#0 atom 1 B`: sdfrust=`2` reference=`1`

### ogb.atom.num_h
- `methyl_radical__h.sdf#0 atom 0 C`: sdfrust=`1` reference=`0`
- `methyl_radical__noh.sdf#0 atom 0 C`: sdfrust=`4` reference=`3`

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
- `first_5K__noh.sdf#8 bond 1`: sdfrust=`0` reference=`5`
- `first_5K__noh.sdf#8 bond 5`: sdfrust=`0` reference=`5`
- `first_5K__noh.sdf#22 bond 7`: sdfrust=`0` reference=`5`
- `first_5K__noh.sdf#29 bond 1`: sdfrust=`0` reference=`5`

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
