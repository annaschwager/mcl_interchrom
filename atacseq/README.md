## ATAC-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the ATAC-seq data from the article X. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to XXX or downloaded from [EGA](https://ega-archive.org/) (EGAD00001002902, EGAD00001002918).

### File pre-processing
Raw sequencing data was processed using the [nf-core/atacseq](https://nf-co.re/atacseq/2.0) pipeline (v2.0) with default parameters. 
Briefly, the files were trimmed to get rid of the sequencing adapters using TrimGalore (v 0.6.7) and aligned to the reference human genome (GRCh38) using BWA (v 0.7.17-r1188). The PCR duplicates were marked using Picard MarkDuplicates function (v 2.27.4). The alignments were filtered for duplicates, unmapped reads, reads not marked as primary alignments, multi mapping reads, reads mapping to the GRCh38 blacklist regions and mitochondrial DNA using SAMtools (v 1.16.1); for soft-clipped reads, reads with more than 4 mismatches and reads with the insert size of >2kb using BAMtools (v 2.5.2). The bigWig files with the ATAC-seq signal scaled to 1 million mapped reads were generated using BEDtools (v 2.30.0) and UCSC BedGraphToBigWig (v 377) and visualized in the IGV Genome Browser. Peak calling was performed using MACS2 (v 2.2.7.1) with the default pipeline parameters.

### Scripts in this folder

**1. diffbind_patients_atac.R** \
This script takes the bam files and the ATAC peaks detected by [MACS2](https://pypi.org/project/MACS2/) in MCL patients samples and naive B cell samples from healthy donors and identifies the differentially acessible sites between these conditions using [Diffbind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html). It than performs the annotation and functional enrichment on the differential peaks.

*Used for:* Figure 2 a b c

**2. diffbind_abemin_atac.R** \
This script takes the bam files and the ATAC peaks detected by [MACS2](https://pypi.org/project/MACS2/) in MCL (GRANTA) and control lymphoblastoid (BLAS) cell lines treated or not with Minnelide (50 nM, 3days) or Abemaciclib (500nM, 7 days) and identifies the differentially acessible sites between these conditions using [Diffbind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html). It than performs the functional annotation of the differentially accessible sites, compares these sites between conditions, and plots the distributions of these sites per chromosome. 

*Used for:* Figure 6 
