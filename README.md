## GWAS-Workflow

Performing a Genome-Wide Association Studies (GWAS) can be a laborious task.
In this reporsitory, we describe the pipeline that we have adapted to carry on a GWAS analysis using **ddRAD-seq-derived SNPs**.
Underneath this README, we illustrate the different steps needed to run the analysis and we provide the codes necessary to reproduce this work.


<br />

<p align="center">
  <img width="500" height="500" src="https://github.com/najlaksouri/GWAS-Workflow/blob/main/Images/Pipeline.png">
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

## Content
1. ddRAD-Sequencing
2. [Data processing](https://github.com/najlaksouri/GWAS-Workflow/tree/main/02.%20Data%20Processing)
    - De-multiplexing: ------------------------> [stacks.sh](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/Stacks.sh)
    - QC filtering: ------------------------> [stacks.sh](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/Trimmomatic.sh)
    - Read mapping: ------------------------> [bwa.mk](https://github.com/najlaksouri/GWAS-Workflow/blob/main/02.%20Data%20Processing/bwa.mk)
3. SNP calling and filtering
    - BCFtools calling
    - GATK calling
    - Freebayes calling
    - Intersection
    - Imputation 
4. GWAS analysis
   - Phenotypic data assessement
   - Population Structure
   - Relatedness
   - Statistical model assessment
   

## Reference


## That's it.... Good luck!
