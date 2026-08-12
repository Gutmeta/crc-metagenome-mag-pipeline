# Species-Level Disease Association Overlap Report

## Definition
Disease-associated species are the top 300 eligible disease-level ranked MetaPhlAn4 species for CRA, CRC and IBD; IBS contributes all eligible ranked species if fewer than 300 are available.
Permutation tests used 10000 fixed-size draws from disease-specific detectable species backgrounds.
IBS association ranking was augmented with external differential species tables from Goll_2020_FMT_IBS, PRJEB34103_IBS.

## Associated Set Sizes
- CRA: 300
- CRC: 300
- IBD: 300
- IBS: 238

## Primary Top300 Overlap
- CRA|CRC: shared=67, expected=49.56, obs/exp=1.35, perm_q=0.006874, significant=True
- CRA|IBD: shared=62, expected=49.98, obs/exp=1.24, perm_q=0.05371, significant=False
- CRA|IBS: shared=35, expected=42.32, obs/exp=0.83, perm_q=1, significant=False
- CRC|IBD: shared=75, expected=48.83, obs/exp=1.54, perm_q=0.0011, significant=True
- CRC|IBS: shared=43, expected=41.53, obs/exp=1.04, perm_q=0.5799, significant=False
- IBD|IBS: shared=61, expected=49.42, obs/exp=1.23, perm_q=0.05371, significant=False

## Multi-Disease Sharing
- CRA|CRC|IBD: shared=19, obs/exp=2.64, perm_q=0.00275
- CRA|CRC|IBS: shared=6, obs/exp=0.98, perm_q=0.7065
- CRA|IBD|IBS: shared=8, obs/exp=1.29, perm_q=0.4384
- CRC|IBD|IBS: shared=15, obs/exp=2.53, perm_q=0.004766
- CRA|CRC|IBD|IBS: shared=0, obs/exp=0.00, perm_q=1

Species shared by >=2 diseases: 247
Species shared by >=3 diseases: 48
Named species shared by >=3 diseases: 30
Species shared by all 4 diseases: 0

## Direction Summary
- CRA|CRC: shared=42, same_direction=34, mixed_direction=8, named=19
- CRC|IBD: shared=41, same_direction=34, mixed_direction=7, named=28
- IBD|IBS: shared=38, same_direction=10, mixed_direction=28, named=13
- CRA|IBD: shared=35, same_direction=11, mixed_direction=24, named=26
- CRC|IBS: shared=22, same_direction=17, mixed_direction=5, named=16
- CRA|IBS: shared=21, same_direction=4, mixed_direction=17, named=13
- CRA|CRC|IBD: shared=19, same_direction=8, mixed_direction=11, named=12
- CRC|IBD|IBS: shared=15, same_direction=5, mixed_direction=10, named=9
- CRA|IBD|IBS: shared=8, same_direction=0, mixed_direction=8, named=4
- CRA|CRC|IBS: shared=6, same_direction=2, mixed_direction=4, named=5

Conclusion: The data support partial species-level overlap, but not the stronger pre-specified broad-overlap criterion across all four diseases.

Important limitation: this analysis supports overlap among disease-associated species, not overlap among strict key species.
