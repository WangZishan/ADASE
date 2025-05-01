# **Allele Specific Expression in Alzheimer's Disease\**


Code for the analysis in this project.<br /><br />
## **./SampleClassification**

**SampleClassification.R**
Classify all samples into Control or AsymAD or AD, based on a rubric considering CERAD, Braak score and Dementia status used in one previous study [https://www.nature.com/articles/s41591-020-0815-6]. 
<br /><br /><br />


## **./ASEVariantIdentification**

Identification of rare NC P/LPs associated with allele specific expression (ASE).


<br /><br /><br />



## **./ASEVariantEnrichmentAtChromosomalBand**

**AR.R**

Plots for the carrier frequency and NC P/LPs count of autosomal recessive (AR) and autosomal dominant (AD) genes across ancestries.

**Distribution_for_ACMG_classification.R**

Plots for the frequency of NC P/LP carriers and count of NC P/LPs across ancestries.

**Distribution_for_genes.R**

Plots for the frequency/count of NC P/LP carriers in each ancestry among the ACMG 59 genes and the top 10% genes (ranked by sums of all defined ancestry frequencies, excluding Mix and Other).

<br /><br /><br />




## **./ASEVariantAssociationWithClinical**

**VariantImpactOnExp.R**

Identification of genes whose expression is affected by related NC P/LPs.

**VariantImpactOnExp_Plot.R**

Volcano plots for genes whose expression is affected by related NC P/LPs.

**PlotPercentileExp.R**

Distribution of percentile expression in a specific cancer at NC P/LP carriers of genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants. Color of node represents variant type. Color of node edge represent ASE enrichment status.

**PlotPercentileExp_DiffExpSplitCount.R**

Count/Proportion of sample-variants across expression splits vs predicted variant function/ASE status for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants.

**PlotPercentileExp_DiffExpSplitCount_GeneInfo.R**

Detailed information for sample-variants of for genes whose expression is significantly/suggestively impacted by NC P/LPs or enriched with significant ASE variants.

<br /><br /><br />




## **./ADAssociatedASEVariantIdentification**

Command of bcftools to extract information of variant of interests from gnomad dataset.

**count.R**

Variant count of predisposing variants in the matched gnomAD ancestry (European of gnomAD is the union of FIN and NFE populations). TCGA population-specific NC P/LPs, exclusively found in a specific TCGA ancestry, are shown as a triangle. Top NC P/LP or top TCGA ancestry-specific NC P/LP, ranked by allele counts in TCGA or gnomAD, was labelled.

**statistic.R**

(Significance of) Correlations of variant frequencies in the matched ancestries between TCGA and gnomAD.


## Citation
If you used or adapted our code in your study, please cite our paper [Under submission].
<br /><br />
## Contact
If you run into issues, you can contact authors: zishan.wang{at}mssm.edu or kuan-lin.huang{at}mssm.edu.
