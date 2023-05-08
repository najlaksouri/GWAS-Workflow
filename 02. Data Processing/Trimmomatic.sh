#!/bin/bash


# define the variables and input files

TRIMMOMATIC_HOME=/home/sie/softwares/Trimmomatic-0.39

############################################################################################################################
				                                                Get started with pool1    
                                                        
                                                        
#POOL1_HOME=/ddRAD-seq/03.Stacks_demultiplexing/pool1
#POOL1_OUT_P=/ddRAD-seq/04.Trimmed_out/pool1/paired
#POOL1_OUT_U=/ddRAD-seq/04.Trimmed_out/pool1/unpaired
############################################################################################################################



mkdir -p ${POOL1_OUT_P} ${POOL1_OUT_U}
for i in {08..29} # for each sample

do
    java -jar ${TRIMMOMATIC_HOME}/trimmomatic-0.39.jar PE -phred33 -threads 2  ${POOL1_HOME}/sample_${i}.1.fq.gz ${POOL1_HOME}/sample_${i}.2.fq.gz \
    ${POOL1_OUT_P}/sample_${i}_R1_paired.fq.gz ${POOL1_OUT_U}/sample_${i}_R1_unpaired.fastq.gz \
    ${POOL1_OUT_P}/sample_${i}_R2_paired.fq.gz ${POOL1_OUT_U}/sample_${i}_R2_unpaired.fastq.gz HEADCROP:7 MINLEN:36

done
