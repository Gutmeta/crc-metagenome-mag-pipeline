# CRA/CRC/IBD/IBS Key Species Overlap V2 Report

## Included disease-arms
- CRA: 6 (Feng_2015_CRC_Austria, Gao_2021_CRC, Lee_2023_CRC, Thomas_2019_CRC_Italy, Yachida_2019_CRC, Zeller_2014_CRC_France)
- CRC: 9 (Feng_2015_CRC_Austria, Gao_2021_CRC, Liu_2022_CRC_China, Thomas_2019_CRC_Italy, Vogtmann_2016_CRC_USA, Wirbel_2019_CRC_Germany, Yachida_2019_CRC, Yu_2017_CRC_China, Zeller_2014_CRC_France)
- IBD: 9 (Braun_2024_Crohn, Bushman_2020_pediatric_Cdiff_IBD, Damman_2015_FMT_UC, Douglas_2018_child_Crohn, Franzosa_2018_IBD, Hall_2017_IBD, He_2017_Crohn, Kumbhari_2024_IBD, Lloyd-Price_2019_HMP2IBD)
- IBS: 1 (Mars_2020_IBS)

## Key species counts
- CRA: Top50=50; strict=0; ranked=1625
- CRC: Top50=50; strict=27; ranked=1650
- IBD: Top50=50; strict=70; ranked=1357
- IBS: Top50=50; strict=30; ranked=465

Overlap tests used 10000 permutations for the main Top50 analysis.
- CRA|CRC: shared=2, obs/exp=1.47, overlap_coeff=0.04, perm_q=1, pass=False
- CRA|IBD: shared=2, obs/exp=1.43, overlap_coeff=0.04, perm_q=1, pass=False
- CRA|IBS: shared=1, obs/exp=0.66, overlap_coeff=0.02, perm_q=1, pass=False
- CRC|IBD: shared=2, obs/exp=1.48, overlap_coeff=0.04, perm_q=1, pass=False
- CRC|IBS: shared=2, obs/exp=1.33, overlap_coeff=0.04, perm_q=1, pass=False
- IBD|IBS: shared=2, obs/exp=1.10, overlap_coeff=0.04, perm_q=1, pass=False

Shared Top50 species in >=2 diseases: 11
Shared Top50 species in >=3 diseases: 0
Shared Top50 species in all 4 diseases: 0

Conclusion: the pre-specified full four-disease claim is not fully supported; report the strongest passing disease-pair and multi-way results instead.
Failed criteria: pairwise_pass_counts={'CRA': 0, 'CRC': 0, 'IBD': 0, 'IBS': 0}; no_passing_3_or_4_way_overlap_containing_CRA_and_IBS; multi_disease_shared_species_count_insufficient(shared3=0,shared4=0)
