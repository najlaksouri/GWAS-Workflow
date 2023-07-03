## Introduction

<div align="justify">


SNPs calling process was carried out using three variant callers:
  - BCFtools v1.7 (Li, 2011) 
  - GATK-HaplotypeCaller v4.2.3.0 (Van der Auwera et al., 2013)
  - Freebayes v1.0.0 (Garrison and Marth, 2012)

Variant calling (SNPs and indels) was conducted in a single-sample mode. Therefore, raw SNPs underwent standard quality filtering based on:
-  **mapping quality (MQ > 40)**
-  **variant quality (QUAL < 30)**
- **depth of eads (DP ≥ 5)**, to remove artifactual calls. 

Clean SNPs from each calling method were merged by position and by reference/alternative alleles into multi-sample VCF files. Then, they were filtered by **call rate > 80%** and residual missing genotypes were imputed with beagle’s default settings. SNPs with minor allele frequency (MAF > 0.05) were selected as a final call set to determine the population structure and marker-trait associations.

</div>


## Citations:

- Garrison E., Marth G. Haplotype-based variant detection from short-read sequencing (2012). arXiv preprint. https://doi.org/10.48550/arXiv.1207.3907

- Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. (2011) Bioinformatics 1;27(21):2987-93. https://doi.org/10.1093/bioinformatics/btr509

- Van der Auwera GA., Carneiro M., Hartl C., Poplin R., Del Angel G., Levy-Moonshine A., Jordan T., Shakir K., Roazen D., Thibault J., et al. From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline (2013). Curr Protoc Bioinformatics. https://doi.org/10.1002/0471250953.bi1110s43
