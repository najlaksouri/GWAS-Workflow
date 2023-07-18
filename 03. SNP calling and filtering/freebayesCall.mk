#####################################################################################
#
# This makefile is to perform raw varinat calling using Freebayes, then filter and
# merge results from different samples
#
# Author: Najla Ksouri
#
#####################################################################################


MAKEFILE=~/ddRADseq/scripts/freebayesCall.mk
MAKE=make -s -f ${MAKEFILE}


# 1. Declare the variable

VAR_DIR=~/ddRAD-seq/08.Variant_calling
FREEBAYES_DIR=${VAR_DIR}/freebayes_call
RAW_VAR=${VAR_DIR}/freebayes_call/raw_call
FILTERED_VAR=${VAR_DIR}/freebayes_call/filtered_call
GENOME_FILE=~/ddRAD-seq/06.Mapping_BWA/ref_genome/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna
INPUT_BAMS=~/ddRAD-seq/07.Post_alignment/srt_coord

samples= 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25R 26 27 28 29 


create_dir:
	echo "create directory"
	mkdir -p ${FREEBAYES_DIR}
	mkdir -p ${RAW_VAR}
	mkdir -p ${FILTERED_VAR}



#################################################################################
#				  reminder freebayes				#
# 										#  
# -min-mapping-quality 40, exlude alignment from analysis if they have a 	#
# mapping quality less than Q, in my case min mapping quality is set at 40	#
#										#
#################################################################################


freebayes:$(samples)
$(samples):
	@echo
	echo "call variants for sample_$@"
	freebayes --min-mapping-quality 40 -f ${GENOME_FILE} ${INPUT_BAMS}/sample_$@.fix_srt.bam > ${RAW_VAR}/raw_free_sample_$@.vcf


freebayes_filter:$(samples)
$(samples):
	echo "filter rawvariants of sample_$@"
	bcftools filter -i "%QUAL > 30 && DP >= 5" ${RAW_VAR}/raw_free_sample_$@.vcf > ${FILTERED_VAR}/filtered_sample_$@.vcf
	echo "compress and index sample_$@"
	bgzip ${FILTERED_VAR}/filtered_sample_$@.vcf
	tabix ${FILTERED_VAR}/filtered_sample_$@.vcf.gz



freebayes_merge:
	mkdir -p ${FREEBAYES_DIR}/merged_call
	echo "merge filtered vcf files"
	bcftools merge --threads 6 ${FILTERED_VAR}/*.vcf.gz -Oz -o ${FREEBAYES_DIR}/merged_call/merged_filtered.vcf.gz



#################################################################################
#				reminder          				#
# 										#		
#  After merging the individual vcf files into a multisample vcf. I have faced  #
# the problem of having a lot of missing called gentoyped, but it was           #
# impossible to know if these missing data are actually missing calls or 	#
# ref homozygote genotype.							#
#										#
# to fix this issue i have used the java package:				#
# http://lindenb.github.io/jvarkit/FixVcfMissingGenotypes.html			#
# this program read a VCF, look back at some BAMS to tells if the missing 	#
# genotypes were homozygotes-ref or not-called. 				#
# If the number of reads is greater than min.depth, then a missing genotype is 	#
# said hom-ref. Herein min.depth was set to 5 					#
# -B: imput.list is txt file containing paths to the bam files			#
# -f Update all fields like DP even if the Genotype is called			#
#################################################################################

FixVcfMissingGenotypes_dir=~/software/jvarkit/dist/

freebayes_fix:
	echo "here we will be using the FixVcfMissingGenotypes java package"
	echo "create the imput list of the bam files"
	find ${INPUT_BAMS} -name "*fix_srt.bam" > ${FREEBAYES_DIR}/merged_call/input.list
	echo "remove FORMAT/GT info because it gives error when fixing missing gentotypes"
	bcftools annotate -x 'FORMAT/GL' ${FREEBAYES_DIR}/merged_call/merged_filtered182.vcf.gz > ${FREEBAYES_DIR}/merged_call/merged_filtered182_mod_format.vcf
	bgzip ${FREEBAYES_DIR}/merged_call/merged_filtered_modified_format.vcf
	@echo "run fix"
	java -jar ${FixVcfMissingGenotypes_dir}/fixvcfmissinggenotypes.jar -f -d 5 -B ${FREEBAYES_DIR}/merged_call/input.list \
	${FREEBAYES_DIR}/merged_call/merged_filtered_modified_format.vcf.gz > ${FREEBAYES_DIR}/merged_call/fixed.merged_filtered.vcf


GATK_HOME=~/software/gatk-4.2.3.0

separate_var:
	echo "keep only snps"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar SelectVariants -R ${GENOME_FILE} -V  ${FREEBAYES_DIR}/merged_call/fixed.merged_filtered.vcf \
  	--select-type-to-include SNP -O ${FREEBAYES_DIR}/merged_call/clean_SNPs.vcf.gz
	bcftools stats ${FREEBAYES_DIR}/merged_call/clean_SNPs.vcf.gz | head -30
	@echo
	@echo
	echo "keep only indels"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar SelectVariants -R ${GENOME_FILE} \
	-V  ${FREEBAYES_DIR}/merged_call/fixed.merged_filtered.vcf --select-type-to-include INDEL -O ${FREEBAYES_DIR}/merged_call/clean_Indels.vcf.gz
	bcftools stats ${FREEBAYES_DIR}/merged_call/clean_Indels.vcf.gz | head -30
