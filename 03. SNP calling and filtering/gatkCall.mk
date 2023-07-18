#####################################################################################
# This makefile is to call variants using GATK-HAPLOTYPECALLER, then filter and
# merge results from different samples
# 
# Author: Najla Ksouri
#
# https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4
#####################################################################################


MAKEFILE=/home/sie/ddRADseq/scripts/gatkCall.mk
MAKE=make -s -f ${MAKEFILE}

PICARD_HOME=~/softwares/
GATK_HOME=~/softwares/gatk-4.2.3.0

#---------------------------------------- 1. Index the genome--------------------------------

#################################################################################	
#			remainder, index 					#
#										#
#										#
# To index the genome, we will need to files (.fai) obtained			#
# from "samtools faidx" and (.dict) form "picard, CreateSequenceDisctionnary"	#
#										#
#################################################################################


index_ref:
	@echo
	echo "create the fai and dict genome files"
	samtools faidx ${REF.FASTA}
	java -jar ${PICARD_HOME}/picard.jar CreateSequenceDictionary -R ${REF.FASTA}


#---------------------------------------- 2. Index BAM files --------------------------------


index_bams:
	@echo
	for i in ${samples};\
	do \
		echo "index the bam alignement of sample_$${i}";\
		samtools index ${INPUT_BAMS}/sample_$${i}.fix_srt.bam;\
	done


#---------------------------------------- 3. run haplotype in vcf mode (variant only) -----------


VCF_MODE=~/ddRAD-seq/08.Variant_calling/gatk_call/vcf_mode
RAW_VCF_CALL=${VCF_MODE}/raw_call
FILTERED_VCF_CALL=${VCF_MODE}/filtered_call


gatkcall:$(samples)
$(samples):
	mkdir -p ${VCF_MODE}
	mkdir -p ${RAW_VCF_CALL}
	mkdir -p ${FILTERED_VCF_CALL}
	echo "call sample_$@"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar HaplotypeCaller -I ${INPUT_BAMS}/sample_$@.fix_srt.bam -R ${REF.FASTA} -O ${RAW_VCF_CALL}/sample_$@.vcf.gz
	echo "filter sample_$@"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar VariantFiltration -V ${RAW_VCF_CALL}/sample_$@.vcf.gz \
	--filter-expression "QUAL < 30.0" --filter-name "filter.QUAL" --filter-expression "MQ < 40.0" --filter-name "map.40" \
	--filter-expression "DP <= 5.0" --filter-name "site_depth" -O ${FILTERED_VCF_CALL}/filtered_sample_$@.vcf.gz
	echo "remove filtered variants in sample_$@"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar SelectVariants -R ${REF.FASTA}  -V ${FILTERED_VCF_CALL}/filtered_sample_$@.vcf.gz --exclude-filtered true \
	-O ${FILTERED_VCF_CALL}/retained_sample_$@.vcf.gz



MERGED_VCF_CALL=${VCF_MODE}/merged_call

merging:
	mkdir -p {MERGED_VCF_CALL}
	bcftools merge --threads 4  ${FILTERED_VCF_CALL}/retained_sample_*.vcf.gz -Oz -o ${MERGED_VCF_CALL}/merged_filtered_Var.vcf.gz
	tabix ${MERGED_VCF_CALL}/merged_filtered_Var.vcf.gz
	@echo
	echo "separate  snps"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar SelectVariants -R ${REF.FASTA} \
	-V ${MERGED_VCF_CALL}/merged_filtered_Var.vcf.gz --select-type-to-include SNP -O ${MERGED_VCF_CALL}/snps_Var.vcf.gz
	@echo 
	echo "separate  indels"
	java -jar ${GATK_HOME}/gatk-package-4.2.3.0-local.jar SelectVariants -R ${REF.FASTA} \
	-V ${MERGED_VCF_CALL}/merged_filtered_Var.vcf.gz --select-type-to-include INDEL -O ${MERGED_VCF_CALL}/indels_Var.vcf.gz
	@echo
	echo "done"


FixVcfMissingGenotypes_dir=~/softwares/jvarkit/dist/

snps_fix:
	#echo "here we will be using the FixVcfMissingGenotypes java package"
	@echo "run fix"
	java -jar ${FixVcfMissingGenotypes_dir}/fixvcfmissinggenotypes.jar -f -d 5 -B ${VCF_MODE}/input.list \
	${MERGED_VCF_CALL}/snps_Var.vcf.gz > ${MERGED_VCF_CALL}/fixed_snps_Var.vcf
	bgzip ${MERGED_VCF_CALL}/fixed_snps_Var.vcf



