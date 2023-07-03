#####################################################################################
# This makefile is to perform raw varinat calling using bcftools, then filter and
# merge results from different samples
#
# Author: Najla Ksouri
#####################################################################################

# 1. Declare the variable

VAR_DIR=/media/sie/TOSHIBA_EXT/Najla/ddRAD-seq/08.Variant_calling
BCFTOOLS_DIR=${VAR_DIR}/bcftools_call
#RAW_VAR=${VAR_DIR}/bcftools_call/raw_call
#FILTRED_VAR=${VAR_DIR}/bcftools_call/filtred_call

GENOME_FILE=~/ddRAD-seq/06.Mapping_BWA/ref_genome/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna
GENOMEFAI=~/ddRAD-seq/06.Mapping_BWA/ref_genome/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fna.fai

INPUT_BAMS=~/ddRAD-seq/07.Post_alignment/srt_coord


#---------------------------------------- 1. Create the index genome--------------------------------


##################################################################################
#			                    reminder samtools, faidx									
#																				                      
# this will create a fasta index with .fai suffix							
#																		
##################################################################################

index:
	echo "create bcftools directory"
	mkdir -p ${BCFTOOLS_DIR}
	echo "index the genome using samtools"   # this will create a file.fai
	samtools faidx ${GENOME_FILE}
	echo ${GENOMEFAI}
	echo "done"

 
 #-----------------------------------  2. bcftools mpileup and call --------------------------------


#################################################################################
#				          reminder bcftools, call 								                      #
#																				                                        #
# mpileup command generate a bcf file containing the gentoype likelihood,		    #
# (count the read coverage and convert the bam to genomic positions)  			    #	
# 																				                                      #
# -Ou: generate an uncompressed output file                                     #
# -f refers to the reference													                          #
# -m: Alternative model for multiallelic and rare-variant calling designed to 	#
#  overcome known limitations in -c calling (the previous version)				      #
# -v flag: output potential variant sites only									                #
# -O z: Output type. “z” means the VCF file is compressed.						          #
#################################################################################


samples=#08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25R 26 27 28 29 


create_dir: 
	mkdir -p ${BCFTOOLS_DIR} 
	mkdir -p ${RAW_VAR} 
	mkdir -p ${FILTRED_VAR}


bcfvariant:$(samples)
$(samples):
	echo "Run raw bcf call for sample_$@" 
	bcftools mpileup -Ou -f ${GENOME_FILE} ${INPUT_BAMS}/sample_$@.fix_srt.bam | bcftools call -mO z  -o ${BCF_RAW_VAR}/rawCall.sample_$@.vcf.gz


 # -----------------------------------  3. bcftools filter and merge ------------------------------


# filtering data based on the 


#################################################################################
#				              reminder bcftools, filter						                  		#
#																				                                        #
# we are goint to filter based on Quality(QUAL), raw read depth (DP) and        #
# average quality mapping (MQ)  								                                #
#                                                                               #
#################################################################################


BCF_RAW_VAR=~/ddRAD-seq/08.Variant_calling/bcftools_call/raw_call
BCF_FILTERED_VAR=~/ddRAD-seq/08.Variant_calling/bcftools_call/filtred_call
BCF_MERGED_VAR=~/ddRAD-seq/08.Variant_calling/bcftools_call/merged_call
MIN-AC-CALL=~/ddRAD-seq/08.Variant_calling/bcftools_call/min-ac

bcf_filter:$(samples)
$(samples):
	mkdir -p ${MIN-AC-CALL}
	echo "filter rawvariants of sample_$@"
	bcftools filter -i "%QUAL > 30 && DP >= 5 && MQ > 40" ${BCF_RAW_VAR}/rawCall.sample_$@.vcf.gz > ${BCF_FILTERED_VAR}/filtered_sample_$@.vcf
	echo "compress and index sample_$@"
	bgzip ${BCF_FILTERED_VAR}/filtered_sample_$@.vcf
	bcftools index ${BCF_FILTERED_VAR}/filtered_sample_$@.vcf.gz
	echo "index raw call for sample$@"
	bcftools index ${BCF_RAW_VAR}/rawCall.sample_$@.vcf.gz
	echo "remove non variant site for sample_$@"
	bcftools view --min-ac=1 ${BCF_FILTERED_VAR}/filtered_sample_$@.vcf.gz -Oz -o ${MIN-AC-CALL}/min_ac_filtered_sample_$@.vcf.gz


 
 merge:
	mkdir -p ${BCF_MERGED_VAR}
	echo "merge all filtered calls"
	bcftools merge --threads 6 ${MIN-AC-CALL}/*vcf.gz -Oz -o ${BCF_MERGED_VAR}/all_merged_filtered.vcf.gz

 
only_var:
	echo "remove non-variant position "
	bcftools view --min-ac=1 ${BCF_MERGED_VAR}/all_merged_filtered.vcf.gz -Oz -o ${BCF_MERGED_VAR}/all_merged_filtered_onlyVar.vcf.gz 


 separate_var:
	echo "keep only snps"
	bcftools filter -i 'TYPE="snp"' ${BCF_MERGED_VAR}/all_merged_filtered_onlyVar.vcf.gz -Oz -o ${BCF_MERGED_VAR}/clean_SNPs.vcf.gz
	@echo
	echo "keep only indels"
	bcftools filter -i 'TYPE="indel"' ${BCF_MERGED_VAR}/all_merged_filtered_onlyVar.vcf.gz  -Oz -o  ${BCF_MERGED_VAR}/clean_Indels.vcf.gz


