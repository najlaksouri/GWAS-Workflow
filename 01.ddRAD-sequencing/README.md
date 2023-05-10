
## Introduction

Double digest RAD-seq (ddRAD-seq) is a flexible and cost-effective strategy thathas emerged as one of the most popular genotyping approaches in plants. It relies on combining two restriction enzymes for library preparation followed by PCR amplification of the template molecules. Briefly, ddRAD-seq is based on breaking the genome into DNA fragments of a certain size using two different restriction enzymes:
  - one common cutter with a short recognition site and
  - a low frequency cutter having a large recognition motif (Aballay et al., 2021). 

This flexible methodology reduces the genotyping complexity and the cost per sample while increasing the number of samples per run. After DNA fragmentation, barcoded adapters are ligated to the end of each DNA strand. Barcoded adaptor is simply a sequence containing the overhang restriction enzyme plus a short oligonucleotide that serves to tag reads from the same sample. Only those falling between both restriction sites and with a specific size range (400-500 bp) are PCR amplified using Illumina indexed primers and subsequently sequenced (Aguirre et al., 2019). An indexed primer is an oligonucleotide with two portions: Illumina primers plus an index (8 bp) which allows the identification of each library (Aguirre et al., 2019).

PstI and MboI restriction enzymes were selected as the best enzyme pair combination
as they produced the highest number of loci with a size range between 400 and 500 bp.
