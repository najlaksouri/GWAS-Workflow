

#####################################################################################
# This makefile is to filter SNPs by call rate and minor allele frequency (MAF), then
# perform the imputation
#
# Author: Najla Ksouri
#####################################################################################



# 1. Declare the variable

GWAS_DIR=~/ddRAD-seq/09.GWAS
MISS_MAF_DIR=${GWAS_DIR}/filter_miss_maf_impute
snps_cv=${MISS_MAF_DIR}/snps_cv.vcf.gz


remove_snps_scaffold:
  echo "create directories"
  mkdir -p ${GWAS_DIR} ${MISS_MAF_DIR}
  echo "remove SNPs in scaffols "
  zless ${snps_cv} | sed '/NW_01802/d' | sed '/NC_014697.1/d' > ${MISS_MAF_DIR}/snps_cv_scaff.vcf
  bgzip ${MISS_MAF_DIR}/snps_cv_scaff.vcf
  tabix -f ${MISS_MAF_DIR}//snps_cv_scaff.vcf.gz


keep_biallelic:
  echo "remove multiallelic SNPs"
  bcftools view -m2 -M2 ${MISS_MAF_DIR}//snps_cv_scaff.vcf.gz -Oz -o ${MISS_MAF_DIR}/biallelic_snps_cv.vcf.gz


missing_call:
  echo "filter by missing call"
  bcftools filter -e 'F_MISSING > 0.2' ${MISS_MAF_DIR}/biallelic_snps_cv.vcf.gz -Oz -o ${MISS_MAF_DIR}/miss_bi_snps_cv.vcf.gz
  echo "snps in Cv:\n" && bcftools stats ${MISS_MAF_DIR}/miss_bi_snps_cv.vcf.gz | head -30

impute:
  echo "impute SNPs"
  beagle gtgl=${MISS_MAF_DIR}/miss_bi_snps_cv.vcf.gz impute=true out=${MISS_MAF_DIR}/bgl_miss_bi_snps_cv

maf:
  echo "filter by maf"
  tabix -f ${MISS_MAF_DIR}/bgl_miss_bi_snps_cv.vcf.gz
  bcftools filter -e 'MAF < 0.05' ${MISS_MAF_DIR}/bgl_miss_bi_snps_cv.vcf.gz -Oz -o ${MISS_MAF_DIR}/maf_bgl_miss_bi_snps_cv.vcf.gz
  echo "snps in Cv:\n" && bcftools stats ${MISS_MAF_DIR}/maf_bgl_miss_bi_snps_cv.vcf.gz | head -30 

