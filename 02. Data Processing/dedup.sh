#!/bin/bash

#========================================================================================================
#    This script is to remove duplicated reads from trimmed reads, before proceeding with the alignment
#	   We use the clone_filter command from Stacks program v2.59
#	   https://catchenlab.life.illinois.edu/stacks/manual/#install
#========================================================================================================

# Define variables
CLONE_FILTER_HOME=/ddRAD-seq/05.Clone_filter
DISCARDED_READS=/ddRAD-seq/05.Clone_filter/discarded_reads
DECLONED_READS=/ddRAD-seq/05.Clone_filter/decloned_reads
POOL1=/ddRAD-seq/04.Trimmed_out/pool1/paired/


# Create directories
mkdir -p ${CLONE_FILTER_HOME}
mkdir -p ${DISCARDED_READS}
mkdir -p ${DECLONED_READS}


#------------------------------STEP 1. Discard PCR duplicated in pool1 ------------------------------#

#############################################
#         reminder, clone_filter          ###
# -1 R1.reads                             ### 
# -2 R2.reads                             ###
# -i format of the input files            ###
# -D capture discarded reads to a file    ###            
#############################################

for i in {08..29}

do
	echo "filtering PCR clones from sample_${i}"
	clone_filter -1 ${POOL1}/sample_${i}_R1_paired.fq.gz -2 ${POOL1}/sample_${i}_R2_paired.fq.gz -i gzfastq -D \
	-o ${DECLONED_READS}/ &> ${DISCARDED_READS}/countclone_sample_${i}.log

done
