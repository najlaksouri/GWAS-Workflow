#!/bin/bash

# Description: This script is used to de-multiplex samples from fastqc files and to remove low quality reads
# In case of paired-end it will generate 4 files: sample1.fastq, sample2,fastq, sample1.rem.fastq and sample1.rem.fastq
# The rem files are those containing low quality reads.
# Author: Najla ksouri nksouri@eead.csic.es

# For more details please check the reference: https://catchenlab.life.illinois.edu/stacks/manual/#procrad

#Declare the variables:

OUTDIR_POOL1=~/ddRAD-seq/03.Stacks_demultiplexing/pool1/
POOL1_R1=~/ddRAD-seq/00.Raw_reads/sar08a29-2020Mel_P1_p1_S1_R1_001.fastq.gz
POOL1_R2=~/ddRAD-seq/00.Raw_reads/sar08a29-2020Mel_P1_p1_S1_R2_001.fastq.gz
barcode_pool1=~/ddRAD-seq/03.Stacks_demultiplexing/Barcodes/barcode_pool1.txt

##################################################################################################################################
# The enzyme pair PstI/MboI was selected as it generated the highest number of loci in concordance with the in silico analysis
# the enzymes used were pstI and Mbol1 however the restriction enzyme available in stacks does not include Mbol1,
# so we replace it by "sau3AI" as they have  the same sequence which is GATC/CTAG
# —1 first input file in a set of paired-end sequences.
# —2 second input file in a set of paired-end sequences.
# c,--clean — clean data, remove any read with an uncalled base.
# q,--quality — discard reads with low quality scores.
# r,--rescue — rescue barcodes and RAD-Tags.
###################################################################################################################################



mkdir -p ${OUTDIR_POOL1} 
process_radtags -E phred33 -1 ${POOL1_R1} -2 ${POOL1_R2} -b ${barcode_pool1} --renz_1 pstI --renz_2 sau3AI -r -q -c -i gzfastq --inline_inline -o ${OUTDIR_POOL1}
