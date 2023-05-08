SNPs calling process was carried out using three variant callers:
  - BCFtools
  - Gatk-HaplotypeCaller
  - freebayes

After performing the raw calling in a single-sample mode, DNA variants were filtered based on 
  - variant quality
  - depth of read
  - mapping quality

Subsequently clean VCF files resulting from each calling pipeline were merged together resulting in a multi-samp
