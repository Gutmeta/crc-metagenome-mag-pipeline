# Top-Journal Robustness Summary

## Prediction performance
- CRA: median ROC-AUC=0.574, mean=0.556, range=0.475-0.615, n=6 cohorts.
- CRC: median ROC-AUC=0.763, mean=0.757, range=0.592-0.837, n=9 cohorts.
- IBD: median ROC-AUC=0.790, mean=0.774, range=0.533-0.929, n=9 cohorts.
- CRA vs CRC: mean delta=0.201, exact permutation q=0.0015.
- CRA vs IBD: mean delta=0.219, exact permutation q=0.0015.
- CRC vs IBD: mean delta=0.018, exact permutation q=0.721.

## Overlap robustness
- CRC|IBD Top300 overlap is supported: shared=75, obs/exp=1.54, q=0.0011.
- CRC|IBD|IBS Top300 three-way overlap is supported: shared=15, obs/exp=2.53, q=0.00477.
- CRC|IBS pairwise is not supported as enriched: shared=43, obs/exp=1.04, q=0.58.
- IBD|IBS pairwise is borderline/not FDR-significant: shared=61, obs/exp=1.23, q=0.0537.
- Four-disease shared species are not observed in the Top300 analysis: shared=0.

## Covariate adjustment
- Age+sex adjusted sensitivity was attempted for 19 evaluable disease arms.
- Age+sex adjusted overlap generated 44 tests; 15 passed the permutation q<0.05 and q95 rule.

Interpretation: the supported claim is narrower than broad four-disease overlap. CRA prediction is weaker than CRC/IBD, while CRC and IBD show moderate microbial-panel discrimination. Clinical-marker equivalence is not asserted from these data.
