# CRC/CRA-Merged AUC and Overlap Figure

Panel A keeps the locked CRA, CRC and IBD prediction panels unchanged, but labels the locked sizes as final K and displays cohort-level AUC comparison q values.
Final panel sizes are CRA k=143, CRC k=18, IBD k=185; median AUCs are CRA=0.574, CRC=0.763, IBD=0.790.
AUC comparison q values are BH-adjusted exact permutation q values from the topjournal robustness run.
- CRA_vs_CRC: mean AUC difference=0.201, exact_permutation_q=0.001499.
- CRA_vs_IBD: mean AUC difference=0.219, exact_permutation_q=0.001499.
- CRC_vs_IBD: mean AUC difference=0.018, exact_permutation_q=0.7212.

Panels B-C merge CRA and CRC disease-associated Top300 species into one CRC/CRA state using the union of the original CRA and CRC associated sets.
Set sizes: CRC/CRA=533, IBD=300, IBS=238.

## Pairwise overlap
- CRC/CRA|IBD: shared=118, obs/exp=1.43, permutation_q=0.0004, significant=True.
- CRC/CRA|IBS: shared=72, obs/exp=1.06, permutation_q=0.2749, significant=False.
- IBD|IBS: shared=61, obs/exp=1.23, permutation_q=0.03453, significant=True.

## Three-way overlap
- CRC/CRA|IBD|IBS: shared=23, obs/exp=2.38, permutation_q=0.0004, significant=True.

Panel C now displays the three pairwise comparisons plus the three-disease intersection together as observed-vs-expected overlap counts.

Important: the overlap panels use broader disease-associated ranked species sets, not the locked prediction panels from Panel A.
