## RNA-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the RNA-seq data from the article X. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to XXX or downloaded from [EGA](https://ega-archive.org/) (EGAD00001002315, EGAD00001002336, EGAD00001002452).

### File pre-processing
Raw sequencing data was processed using the [nf-core/rnaseq pipeline (v3.10.1)](https://nf-co.re/rnaseq/3.10.1) with default parameters unless stated otherwise. Briefly, the files were trimmed from the sequencing adapters using TrimGalore (v 0.6.7) and aligned to the reference human genome (GRCh38) using STAR (v 2.6.1d). The reads mapping to different genomic features were then quantified by Salmon (v 1.9.0).

### Scripts in this folder 
**1. diff_expression_MCL.R** \
This script takes the count tables obtained at the pre-processing step, performs quality checks, low counts filtering and differential expression analysis with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). It than generates the plots of the differentially expressed genes mapped to their chromosomal locations and the barplots with the numbers of DEGs per chromosome in different comparisons (raw and normalised to the chromosme size or the gene number per chromosome taken from the *chromSizes.csv*). Finally, it performs the Gene Ontology enrichment analysis for the deregulated genes on chromosome 19 and exports the expression data for the ABC analysis (see *chipseq* folder).\

*Data analysed:* B cells from MCL patients (5 sampels from EGAXXX, 4 samples sequenced for this study and deposited to EGAXXX), control naive B cells (6 samples from EGAD00001002315), GRANTA-519 MCL cells (3 samples, GEOXXX).

*Used for:* Figure 3 a,b,c 

**2. diff_expression_abe_min_cells.R** \
This script takes the count tables obtained at the pre-processing step, performs quality checks, low counts filtering and differential expression analysis with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). It than generates the volcano plots with the differentially expressed genes and generates the barplots with the numbers of DEGs per chromosome in different comparisons (raw and normalised to the chromosme size or the gene number per chromosome taken from the *chromSizes.csv*). Finally, it performs the Gene Ontology enrichment analysis for the genes differentially expressed following the treatment.\

*Data analysed:* MCL (GRANTA-519) and control (BLAS) cells, treated with 50nM Minnelide for 3 days or 500mkM Abemaciclib for 7 days with the corresponding non-treated controls (GEOXXX). 

*Used for:* Figure 7 a,b,c,d,e,f

**2. diff_expression_abe_min_patient.R** \
This script performs the same analysis as *diff_expression_abe_min_cells.R* for the cells from an MCL patient.\

*Data analysed:* PBMCs from an MCL patient in a leukimic phase, treated with with 25nM/50nM Minnelide for 3 days or 500mkM Abemaciclib for 7 days with the corresponding non-treated controls (EGAXXX). 

*Used for:* Figure 7 a,b,c,d,e,f
