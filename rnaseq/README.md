## RNA-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the RNA-seq data. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to [GEO](https://www.ncbi.nlm.nih.gov/geo/) (GSE291376), or downloaded from [EGA](https://ega-archive.org/) (EGAD00001002315, EGAD00001002336, EGAD00001002452).

### File pre-processing
Raw sequencing data was processed using the [nf-core/rnaseq pipeline (v3.10.1)](https://nf-co.re/rnaseq/3.10.1) with default parameters unless stated otherwise. Briefly, the files were trimmed from the sequencing adapters using TrimGalore (v 0.6.7) and aligned to the reference human genome (GRCh38) using STAR (v 2.6.1d). The reads mapping to different genomic features were then quantified by Salmon (v 1.9.0).

### Scripts in this folder 
**1. diff_expression_MCL.R** \
Performs DESeq2 comparisons of primary MCL samples and GRANTA-519 cells against naïve B-cell controls, with downstream chromosome enrichment, chr19-focused analyses, GO enrichment, positional GSEA, and visualization of chromosome 19 transcriptional patterns.

*Data analysed:* B cells from MCL patients (5 sampels from EGAD00001002336, 4 samples sequenced for this study), control naive B cells (6 samples from EGAD00001002315), GRANTA-519 MCL cells (3 samples, GSE291376).

**2. diff_expression_abe_min_cells.R** \
Differential expression analysis of MCL cell lines following Minnelide and Abemaciclib treatment. Includes DESeq2 analysis, PCA, gene-set enrichment, treatment-response comparisons, and assessment of transcriptional reversal.

*Data analysed:* MCL (GRANTA-519) and control (BLAS) cells, treated with 50nM Minnelide for 3 days or 500mkM Abemaciclib for 7 days with the corresponding non-treated controls (GSE291376). 


**3. diff_expression_abe_min_patient.R** \
This script performs the same analysis as *diff_expression_abe_min_cells.R* for the cells from an MCL patient.

*Data analysed:* PBMCs from an MCL patient in a leukimic phase, treated with with 25nM/50nM Minnelide for 3 days or 500mkM Abemaciclib for 7 days with the corresponding non-treated controls. 

