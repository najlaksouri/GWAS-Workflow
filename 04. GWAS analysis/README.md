
# Phenotypic traits evaluation

In this study, we used a discovery panel of 90 peach accessions in order to dissect the genetic architecture of 16 fruit-related traits:

<table >
 <tr>
    <td> Harvest date (HvD; Julian days) </td>
    <td> Fruit weight (FW; grams) </td>
    <td> Flesh firmness (FF; newton) </td>
    <td> soluble solids content (SSC; °Brix) </td>
    <tr>
    <td> Titratable acidity (TA; grams malic acid/100 g flesh weight) </td>
    <td>  Ripening index (RI; SSC/TA) </td>
    <td>  Vitamin C (Vit C; mg of ascorbic acid/100 g flesh weight) </td>
    <td>  Total phenolics (Phen; mg of gallic acid equivalents/100 g flesh weight) </td>
    <tr>
    <td>  Contents of flavonoid (Flvs; catechin equivalents/100 g flesh weight) </td>
    <td>  Anthocyanin (ACNs; cyanidin-3-glucoside/kg flesh weight) </td>
    <td>  Sucrose (Suc; g/kg flesh weight) </td>
    <td>  Glucose (Glu; g/kg flesh weight) </td>
    <tr>
    <td>  Fructose (Fruc; g/kg flesh weight) </td>
    <td>  Sorbitol (SRB; g/kg flesh weight) </td>
    <td>  Total sugars (TS; g/kg flesh weight) </td>
    <td> Relative antioxidant capacity (RAC; μg TE/g flesh weight) </td>
    </tr>
 </tr>
</table>
            
Only traits with **H<sup>2</sup>  > 0.5** were considered for association analysis. Consequently, contents of glucose, fructose, sucrose and total sugars were discarded from the subsequent analysis.


<p align="center">
  <img  src="https://github.com/najlaksouri/GWAS-Workflow/blob/main/04.%20GWAS%20analysis/Figure1_600.jpg">
  
</p>



## Marker-trait associations 
For association mapping, seven statistical models, ranging from single to multi-locus, were simultaneously tested in GAPIT v3.1.065:
- General linear model (GLM)
- Mixed linear model (MLM)
- Compressed MLM (CMLM)
- Settlement of MLMs under progressively exclusive relationship (SUPER).
- Multiple loci mixed linear model (MLMM)
- Fixed and random model circulating probability unification (FarmCPU)
- Bayesian-information and linkage-disequilibrium iteratively nested keyway (BLINK)
