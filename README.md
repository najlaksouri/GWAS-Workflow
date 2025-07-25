## GWAS-Workflow

Performing a Genome-Wide Association Studies (GWAS) can be a laborious task.
In this repository, we describe the pipeline that we have adapted to carry on a GWAS analysis using **ddRAD-seq-derived SNPs**.
Underneath this README, we illustrate the different steps needed to run the analysis and we provide the codes necessary to reproduce this work.


<br />

<p align="center">
  <img width="600" height="600" src="./04.%20GWAS%20analysis/Pipeline.png">
  
</p>


## Prerequisites
- [Stacks](https://catchenlab.life.illinois.edu/stacks/)
- [BWA-Mem](https://github.com/lh3/bwa)
- [SAMtools](https://bioinformaticsreview.com/20210404/installing-samtools-on-ubuntu/)
- [BCFtools](https://samtools.github.io/bcftools/)
- [GATK-HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
- [Freebayes](https://github.com/freebayes/freebayes)
- [Plink](https://www.cog-genomics.org/plink/)
- [R program](https://cran.r-project.org/) and [RStudio Desktop](https://posit.co/download/rstudio-desktop/):
  - ggplot
  - GAPIT
  - CMplot
  - SmartPCA

## Content
**Step 1**. [ddRAD-Sequencing](https://github.com/najlaksouri/GWAS-Workflow/tree/main/01.ddRAD-sequencing) 

**Step 2**. [Data processing](https://github.com/najlaksouri/GWAS-Workflow/tree/main/02.%20Data%20Processing) 
- De-multiplexing: ------------------------> [stacks.sh](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/stacks.sh)
- QC filtering: ------------------------> [trimmomatic.sh](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/trimmomatic.sh)
- Remove duplicated reads ------------------------> [dedup.sh](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/dedup.sh)
- Reads mapping: ------------------------> [align_bwa.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/align_bwa.mk) 
    
**Step 3**. [Variant calling and filtering](https://github.com/najlaksouri/GWAS-Workflow/tree/main/03.%20SNP%20calling%20and%20filtering) 
- BCFtools calling: ------------------------> [bcftoolsCall.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/03.%20SNP%20calling%20and%20filtering/bcftoolsCall.mk)
- Freebayes calling: ------------------------> [freebayesCall.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/03.%20SNP%20calling%20and%20filtering/freebayesCall.mk)
- GATK calling: ------------------------> [gatkCall.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/03.%20SNP%20calling%20and%20filtering/gatkCall.mk)
- MAF filtering anf imputation: ------------------------>  [maf_imputation.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/03.%20SNP%20calling%20and%20filtering/maf_imputation.mk)
    
**Step 4**. [GWAS analysis](https://github.com/najlaksouri/GWAS-Workflow/tree/main/04.%20GWAS%20analysis)
   - Phenotypic data assessement ------------------------> [trait transformations](https://najlaksouri.github.io/GWAS-Workflow/04.%20GWAS%20analysis/phenotypic-assessment.html)
   - PCA analysis
   - Relatedness
   - Statistical model assessment
   

## Reference

Ksouri, N., Sánchez, G., Font i Forcada, C. et al. A reproducible ddRAD-seq protocol reveals novel genomic association signatures for fruit-related traits in peach. Plant Methods 21, 101 (2025). https://doi.org/10.1186/s13007-025-01415-3

## That's it.... Good luck!
