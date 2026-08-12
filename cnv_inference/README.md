## RNA-seq-based copy-number analysis

This folder contains the scripts used for copy-number alteration and allelic-imbalance analyses of the RNA-seq data.

### File pre-processing
RNA-seq BAM files were processed with BAFExtract to obtain SNP-level B-allele frequency information for CaSpER analysis. SNPs were extracted using a minimum mapping quality of 50, minimum coverage of 20 reads, minimum alternative allele count of 4, and minimum minor allele frequency of 0.1.

### Scripts in this folder
**1. cnv_prep.sh**
Generates BAFExtract SNP files from RNA-seq BAM files for downstream CaSpER analysis.
\
\
**2. cnv_inference.R**
Infers large-scale copy-number alterations from RNA-seq data using CaSpER, with naïve B cells as controls, and visualizes arm-level CNA profiles. The script also evaluates RNA-derived allelic imbalance across chromosome 19p and around the chr19 probe region using B-allele frequencies.

*Data analysed*: B cells from MCL patients (5 samples from EGAD00001002336 and 4 samples sequenced for this study), control naïve B cells (6 samples from EGAD00001002315), and GRANTA-519 MCL cells (3 samples, GSE291376).
