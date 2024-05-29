#!/bin/bash

[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

### Introducing variables for paths ###
path_to_project="home/anna/data1/NGS_data_current/DAG_MCL/chip/nfcore_h3k27ac"
path_to_bwa_files="$path_to_project""/bwa/mergedLibrary/"
path_to_macs2_files="$path_to_bwa_files""narrowPeak/"
path_to_rose_files="$path_to_bwa_files""rose/"

samples=(MCL1 MCL2 MCL3 MCL4 MCL5 Bnaive1 Bnaive2 Bnaive3 Bgerminal1 Bgerminal2)
mcl=(MCL1 MCL2 MCL3 MCL4 MCL5)
naive=(Bnaive1 Bnaive2 Bnaive3)
germinal=(Bgerminal1 Bgerminal2)

### Calling narrowPeaks with MACS2 ###
for i in ${samples[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}.mLb.clN.sorted.bam \
                   -c "$path_to_bwa_files"${i}_input.mLb.clN.sorted.bam \
                   -n ${i} -g hs \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_macs2.log 
done

### Merging replicates with CHIP-R ###
chipr -i "$path_to_macs2_files"MCL1_peaks.narrowPeak "$path_to_macs2_files"MCL2_peaks.narrowPeak "$path_to_macs2_files"MCL3_peaks.narrowPeak "$path_to_macs2_files"MCL4_peaks.narrowPeak "$path_to_macs2_files"MCL5_peaks.narrowPeak \
      -m 2 -o "$_path_to_chipr"mcl_chipr_all.gff &

chipr -i "$path_to_macs2_files"Bnaive1_peaks.narrowPeak "$path_to_macs2_files"Bnaive2_peaks.narrowPeak "$path_to_macs2_files"Banive3_peaks.narrowPeak \
      -m 2 -o "$_path_to_chipr"bnaive_chipr_all.gff &

chipr -i "$path_to_macs2_files"Bgerminal1_peaks.narrowPeak "$path_to_macs2_files"Bgerminal2_peaks.narrowPeak \
      -m 2 -o "$_path_to_chipr"bgerminal_chipr_all.gff &

### Calling superenhancers with ROSE ###
for i in ${mcl[*]};do
    python ROSE_main.py -g hg38 -i mcl_chipr_all.gff \
                                -r "$path_to_bwa_files"${i}.mLb.clN.sorted.bam \
                                -o "$path_to_rose_files"rose_output_${i} \
                                -c "$path_to_bwa_files"${i}_input.mLb.clN.sorted.bam \
                                -t 2000 -s 12500 &
done

for i in ${naive[*]};do
    python ROSE_main.py -g hg38 -i bnaive_chipr_all.gff \
                                -r "$path_to_bwa_files"${i}.mLb.clN.sorted.bam \
                                -o "$path_to_rose_files"rose_output_${i} \
                                -c "$path_to_bwa_files"${i}_input.mLb.clN.sorted.bam \
                                -t 2000 -s 12500 &
done

for i in ${germinal[*]};do
    python ROSE_main.py -g hg38 -i bgerminal_chipr_all.gff \
                                -r "$path_to_bwa_files"${i}.mLb.clN.sorted.bam \
                                -o "$path_to_rose_files"rose_output_${i} \
                                -c "$path_to_bwa_files"${i}_input.mLb.clN.sorted.bam \
                                -t 2000 -s 12500 &
done


### Software versions ###
# Python 2.7.18
# macs2 2.2.9.1
# chipr 1.1.0
# rose2 1.1.0

