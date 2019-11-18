# Quantifiction of short-fragmented cfDNA to infer circulating tumor DNA and chromosomal instability detection

Tumor-derived cell-free DNA (cfDNA) or circulating tumor DNA (ctDNA) can be detected in the plasma of patients. 
Analysis of ctDNA provide an alternative source of information for non-invasively track tumor progression and response throughout course of therapy. However, larger quantities of cfDNA from noncancerous cells dilute the overall concentration of ctDNA. Recent bioinformatics tools are mostly capable of detecting DNA mutations (SNVs/indels and CNV) from tumor tissue, but are inferior to detect low tumor DNA contribution in cfDNA sample. Developing a tool that considers short DNA fragmentation, which differentiates ctDNA from non-tumor cfDNA, will improve detection sensitivity [1].

We develop here a novel method to evaluate genome-wide short DNA fragmentation from low-coverage whole-genome sequencing data of cfDNA. The fraction of short DNA fragments, having a length of 150 bases or less, is
calculated for each genomic non-overlapping window. A z-score statistic is used to determine if the short DNA fragment in each bin would be significantly higher or lower compared to a group of healthy individuals. Our method reports genomic regions with deviant fragmentation patterns which could be used to infer copy-number alterations in tumor cells. Quantitative measures of fragmentation deviation, referred to as genome-wide z-scores[2]. 

Prerequisite
1. R 3.5.1





## Reference
1. Mouliere, Florent, et al. "Enhanced detection of circulating tumor DNA by fragment size analysis." Science translational medicine 10.466 (2018): eaat4921.
2. Vanderstichele, Adriaan, et al. "Chromosomal instability in cell-free DNA as a highly specific biomarker for detection of ovarian cancer in women with adnexal masses." Clinical Cancer Research 23.9 (2017): 2223-2231.
