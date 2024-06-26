## ChIP-seq analysis 
This folder contains the scripts used for the downstream analysis and visualisation of the ChIP-seq data from the article X. \
The raw sequencing data was generated in the CNRS UMR9018 and deposited to XXX or downloaded from [EGA](https://ega-archive.org/) (EGAD00001001502, EGAD00001001519, EGAD00001002397).

### File pre-processing
Raw sequencing data was processed using the [nf-core/chipseq](https://nf-co.re/chipseq/2.0.0) pipeline (v2.0.0) with default parameters unless mentioned otherwise. 
Briefly, the files were trimmed to get rid of the sequencing adapters using TrimGalore (v 0.6.7) and aligned to the reference human genome (GRCh38) using BWA (v 0.7.17-r1188).
The PCR duplicates were marked using Picard MarkDuplicates function (v 2.27.4). The alignments were filtered for duplicates, unmapped reads, multi mapping reads, reads mapping 
to the GRCh38 blacklist regions and reads mapping for different chromosomes using SAMtools (v 1.15.1); for reads with more than 4 mismatches and the insert size of >2kb using 
BAMtools (v 2.5.2). The bigWig files with the ChIP-seq signal scaled to 1 million mapped reads were generated using BEDtools (v 2.30.0) and UCSC BedGraphToBigWig (v 377).

### Scripts in this folder

**1. diffenrich_h3k27ac.R** \
This script takes the H3K27ac peaks identified by [MACS2](https://pypi.org/project/MACS2/) and identifies the sites differentially enriched for H3K27ac across control and MCL samples using diffbind (insert ref). 
Put here:
diffenrich, matrices for H3K27aC dif enriched peaks genomewide and per chromosome \
*Used for:* Figure 2 a \
Put here or in another script: H3K27ac peaks per chr 

**2. ABC.sh** \
This script uses the activity-by-contact model of enhancer–promoter regulation (ABC) to detect the enhancers and their regulatory elements based on ATAC-seq, H3K27ac ChIP-seq, RNA-seq and HiC data. The source code of the model is taken from [broadinstitute/ABC-Enhancer-Gene-Prediction](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/tree/master). 

**3. rnaseq_with_abc.R** \
This script takes the results from the ABC model, the results from the RNA-seq differential expression analysis (*see rnaseq folder*) and the results from the DiffBind analysis of the H3K27aC ChIP-seq data to plot the deregulated genes predicted to be regulated by enhancer or promotors that are significantly more active (enriched for the H3K27Ac signal) in MCL patients B cells vs control naive B cells from healthy donors. \
*Used for:* Figure 2 f, g 

**4. se.sh** \
This script uses [MACS2](https://pypi.org/project/MACS2/) to detect narrowPeaks from the aligned data for the H3K27Ac histone mark followed by ChIP-R ([rhysnewell/ChIP-R](https://github.com/rhysnewell/ChIP-R)) to detect the consensus peaksets.
The consensus peaks are fed to the rank ordering of super-enhancers (ROSE) algortihm ([stjude/ROSE/](https://github.com/stjude/ROSE/tree/master)) for super-enahcer detection. \
*Used for:* Figure 1 d, e 

**5. se_consensus_intersections.R** \
This script takes the output of the ROSE (from *se.sh*), detects the consensus SEs within each condition and the overlaps between the consensus SEs across conditions. It than associates the SEs to the nearest genes and performs over-representation enrichment analysis on the resulting gene sets. \
*Used for :* Figure where SEs, supplementary with SEs 



