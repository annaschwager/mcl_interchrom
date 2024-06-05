#!/bin/bash

[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/"
path_to_atac="/home/anna/data1/Projects/MCL/atac_dag/results/bwa/mergedReplicate/"
path_to_h3k27ac="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/"
path_to_expression="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/expression/"
path_to_hic="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/average_hic_data/"
path_to_ref="ref/"
samples=(naive_blood germinal)

###STEP 1###
for i in ${samples[*]};do
    macs2 callpeak \
    -t "$path_to_atac"${i}.mRp.clN.sorted.bam \
    -f BAM \
    -n ${i} \
    -g hs \
    -p .1 \
    --call-summits \
    --outdir ABC_output/${i}
done

for i in ${samples[*]};do
    bedtools sort -faidx "$path_to_ref"GRCh38.p13.genome.fa.sizes -i ABC_output/${i}/${i}_peaks.narrowPeak > ABC_output/${i}/${i}_peaks.narrowPeak.sorted
done

for i in ${samples[*]};do
    python src/makeCandidateRegions.py \
      --narrowPeak  ABC_output/${i}/${i}_peaks.narrowPeak.sorted \
      --bam "$path_to_atac"${i}.mRp.clN.sorted.bam \
      --outDir ABC_output/${i} \
      --chrom_sizes "$path_to_ref"GRCh38.p13.genome.fa.sizes \
      --regions_blocklist "$path_to_ref"hg38-blacklist.v2.bed \
      --regions_includelist "$path_to_ref"RefSeqCurated.170308.bed.CollapsedGeneBounds.hg38.TSS500bp.bed \
      --peakExtendFromSummit 250 \
      --nStrongestPeaks 150000 
done

###STEP 2###

files_MCL="DAGfiles/chip/MCL1_EGAF00000887679.mLb.clN.sorted.bam,DAGfiles/chip/MCL2_EGAF00000887736.mLb.clN.sorted.bam,DAGfiles/chip/MCL3_EGAF00000887744.mLb.clN.sorted.bam,DAGfiles/chip/MCL4_EGAF00000887808.mLb.clN.sorted.bam,DAGfiles/chip/MCL5_EGAF00000887809.mLb.clN.sorted.bam"
files_germinal="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/germinal/Bgerminal2_EGAF00000768436.mLb.clN.sorted.bam,/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/germinal/Bgerminal1_EGAF00000768688.mLb.clN.sorted.bam"
files_blood_naive="/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/naive_blood/Bnaive1_EGAF00000769006.mLb.clN.sorted.bam,/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/naive_blood/Bnaive2_EGAF00000769220.mLb.clN.sorted.bam,/home/anna/data1/Projects/MCL/ABCmodel/ABC-Enhancer-Gene-Prediction-master/DAGfiles/chip/naive_blood/Bnaive3_EGAF00000768750.mLb.clN.sorted.bam"

for i in ${samples[*]};do
      python src/run.neighborhoods.py \
       --candidate_enhancer_regions ABC_output//${i}/${i}_peaks.narrowPeak.sorted.candidateRegions.bed \
       --genes "$path_to_ref"gencode.v39.annotation.bed \
       --H3K27ac files_${i}
       --ATAC "$path_to_atac"${i}.mRp.clN.sorted.bam \
       --expression_table "$path_to_expression"${i}_norm_counts_averagevalues.csv \
       --chrom_sizes "$path_to_ref"neighbours.GRCh38.p13.genome.fa.sizes \
       --ubiquitously_expressed_genes "$path_to_ref"UbiquitouslyExpressedGenesHG19.txt \
       --cellType GM12878 \
       --outdir ABC_output/${i}/Neighborhoods/
done

###STEP 3###

for i in ${samples[*]};do
      python src/predict.py \
      --enhancers ABC_output/${i}/Neighborhoods/EnhancerList.txt \
      --genes ABC_output/${i}/Neighborhoods/GeneList.txt \
      --HiCdir "$path_to_hic"average_hic.v2.191020 \
      --chrom_sizes "$path_to_ref"neighbours.GRCh38.p13.genome.fa.sizes \
      --hic_resolution 5000 \
      --scale_hic_using_powerlaw \
      --threshold .02 \
      --cellType GM12878 \
      --outdir ABC_output/${i}/Predictions/ \
      --make_all_putative
done

### Software versions ###
# Python 2.7.18
# macs2 2.2.9.1
# bedtools v2.30.0
