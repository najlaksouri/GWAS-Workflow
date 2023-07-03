
## Introduction 
<div align="justify">

Double digest RAD-seq (ddRAD-seq) is a flexible and cost-effective strategy that has emerged as one of the most popular genotyping approaches in plants. It relies on combining two restriction enzymes for library preparation followed by PCR amplification of the template molecules. 
Briefly, ddRAD-seq is based on breaking the genome into DNA fragments of a certain size using two different restriction enzymes:
  - one common cutter with a short recognition site and
  - a low frequency cutter having a large recognition motif (Aballay et al., 2021).


In our study, *<ins> PstI </ins>* and *<ins> MboI </ins>* restriction enzymes were selected as the best enzyme pair combination as they produced the highest number of loci with a size range between 400 and 500 bp.


This flexible methodology reduces the genotyping complexity and the cost per sample while increasing the number of samples per run. After DNA fragmentation, barcoded adapters are ligated to the end of each DNA strand. Barcoded adaptor is simply a sequence containing the overhang restriction enzyme plus a short oligonucleotide that serves to tag reads from the same sample.

Only those falling between both restriction sites and with a specific size range (300-400 bp) are PCR amplified using Illumina indexed primers and subsequently sequenced (Aguirre et al., 2019). An indexed primer is an oligonucleotide with two portions: Illumina primers plus an index (8 bp) which allows the identification of each library (Aguirre et al., 2019).
Ligated fragments from 24 samples were subsequently pooled together and were PCR amplified with indexed primers to tag each pool.Finally, paired-end reads (250 bp) were generated on an <ins>Illumina NovaSeq 6000</ins> instrument at CIMMYT, Mexico



**Ps**:
Raw sequence reads were submitted to the **European Nucleotide Archive (ENA)** under the project reference **PRJEB62784**.
</div>


## Citations 

- Aballay, M.M., Aguirre, N.C., Filippi, C.V. et al. Fine-tuning the performance of ddRAD-seq in the peach genome. Sci Rep 11, 6298 (2021). https://doi.org/10.1038/s41598-021-85815-0
- Aguirre, N.C.; Filippi, C.V.; Zaina, G.; Rivas, J.G.; Acuña, C.V.; Villalba, P.V.; García, M.N.; González, S.; Rivarola, M.; Martínez, M.C.; et al. Optimizing ddRADseq in Non-Model Species: A Case Study in Eucalyptus dunnii Maiden. Agronomy 2019, 9, 484. https://doi.org/10.3390/agronomy9090484



