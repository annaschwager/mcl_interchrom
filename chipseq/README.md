## ChIP-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the ChIP-seq data from the article X. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to XXX or downloaded from [EGA](https://ega-archive.org/) (EGAD00001001502, EGAD00001001519, EGAD00001002397).

### File pre-processing
Raw sequencing data was processed using the [nf-core/chipseq](https://nf-co.re/chipseq/2.0.0) pipeline (v2.0.0) (Ewels et al. 2022) with default parameters unless mentioned otherwise. 
Briefly, the files were trimmed to get rid of the sequencing adapters using TrimGalore (v 0.6.7) and aligned to the reference human genome (GRCh38) using BWA (v 0.7.17-r1188).
The PCR duplicates were marked using Picard MarkDuplicates function (v 2.27.4). The alignments were filtered for duplicates, unmapped reads, multi mapping reads, reads mapping 
to the GRCh38 blacklist regions and reads mapping for different chromosomes using SAMtools (v 1.15.1); for reads with more than 4 mismatches and the insert size of >2kb using 
BAMtools (v 2.5.2). The bigWig files with the ChIP-seq signal scaled to 1 million mapped reads were generated using BEDtools (v 2.30.0) and UCSC BedGraphToBigWig (v 377).

### Scripts in this folder
**1. macs2_se.sh** \
This script uses MACS2 to detect narrowPeaks from the aligned data for the H3K27Ac histone mark followed by the ROSE algorithm to predict super-enhancers.
