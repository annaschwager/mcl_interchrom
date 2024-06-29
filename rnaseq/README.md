## RNA-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the RNA-seq data from the article X. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to XXX or downloaded from [EGA](https://ega-archive.org/) (EGAD00001002315, EGAD00001002336, EGAD00001002452).

### File pre-processing


### Scripts in this folder 
**1. diff_expression_MCL.R** \
This script takes the count tables obtained at the previous step, performs quality checks, low counts filtering and differential expression analysis with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). It than generates the plots of the differentially expressed genes mapped to their chromosomal locations and the barplots with the numbers of DEGs per chromosome in different comparisons. Finally, it performs the Gene Ontology enrichment analysis for the deregulated genes on chromosome 19 and exports the expression data for the ABC analysis (see *chipseq* folder).

*Used for:* Figure 3 a b 

