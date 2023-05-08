######################################################################
# This makefile is to perform paired-end alignment using 
# the bwa-mem v0.7.17, then convert SAM to BAM and sort by coordinates
#
# Author: Najla Ksouri
######################################################################


MAKEFILE=/ddRADseq/scripts/bwa.mk
MAKE=make -s -f ${MAKEFILE}


#------------------------------ 1. download and index genome ------------------------------

# Download the reference genome from NCBI
MAPPING_HOME=/media/sie/TOSHIBA_EXT/Najla/ddRAD-seq/06.Mapping_BWA
GENOME_HOME=/media/sie/TOSHIBA_EXT/Najla/ddRAD-seq/06.Mapping_BWA/ref_genome/


download:
	mkdir -p ${MAPPING_HOME} ${GENOME_HOME}
	@echo "download the reference genome from NCBI"
	@wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/346/465/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz \
	-P ${GENOME_HOME}



#############################################
#         reminder, wget                  ###
# -P: indicate the output directory       ###
# 					  ###
#############################################


# indexing the reference genome
GENOME_FILE_GZ=/media/sie/TOSHIBA_EXT/Najla/ddRAD-seq/06.Mapping_BWA/ref_genome/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.gz 
GENOME_FILE=/media/sie/TOSHIBA_EXT/Najla/ddRAD-seq/06.Mapping_BWA/ref_genome/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna
index:
	@echo "indexing the genome"
	@echo
	@gunzip ${GENOME_FILE_GZ} 
	@bwa index ${GENOME_FILE}




DECLONED=/ddRAD-seq/05.Clone_filter/decloned_reads


#------------------------------ 2. Mapping reads ------------------------------


## Mapping the reads using bwa-mem algorithm and multiple cores

p1=08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29


#################################################################
#		reminder, bwa-mem			        #
#								#
# -M and -R are needed for GATK compatibility			#
# -M: if a read is splitted (different parts maps to different	#
# regions, marks all parts other than main as secondary		#
# alignment then GATK which ignores secondary alignment).	#
# -R; add a read groupe RG and sample name tag SM.		#
# BWA MEM command to define the Read Group:			#
# -R '@RG\tID:.....\tSM:......\tLB:......\tPL:......'		#
# ID:unique id of a collection of reads, SM:Sample name;	#
# LB:library name; PL: sequencing platform			#
# the most important components are ID,SM AND LB 		#
# Exple: @RG ID:flowcell1.lane1	PL: ILLUMINA  LB:pool1 	SAM:dad #
#								#
#################################################################

.PHONY: reads $(p1)
pool1: $(p1)
$(p1):
	echo "mapping the pair-reads sample_$@_R1 and sample_$@_R2"
	bwa mem -t 8 -M -R "@RG\tID:sar_$@\tSM:sample_$@\tLB:pool1\tPL:ILLUMINA" ${GENOME_FILE} ${DECLONED}/sample_$@_R1_paired.1.fq.gz ${DECLONED}/sample_$@_R2_paired.2.fq.gz \
	2> ${MAPPING_HOME}/bwa.err.pool1 > ${MAPPING_HOME}/sample_$@.sam
