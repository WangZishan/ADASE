# **Allele Specific Expression in Alzheimer's Disease\**
<br /><br />

## **Introduction**
If you used or adapted our code in your study, please cite our paper [Under submission]. 

If you run into issues, you can contact authors: zishan.wang{at}mssm.edu or kuan-lin.huang{at}mssm.edu.

Five folders were listed, with each containing R code for a specific analysis. Each folder includes a corresponding explanatory section as following.
<br /><br /><br />

## **SampleClassification**
**SampleClassification.R**

Classify all samples into Control or AsymAD or AD, based on a rubric considering CERAD, Braak score and Dementia status used in one previous study [https://www.nature.com/articles/s41591-020-0815-6].
<br /><br /><br />


## **ASEVariantIdentification**

**1_BinomialTest.R**

Calculate the P-value assessing whether one variant exhibite allele specific expression based on binomial test.

**2_ASEVariantIdentification.R**

ASE detection from bi-allelic heterozygous variants of autosomal chromosomes that passed stringent quality control in each sample.

<br /><br /><br />



## **ASEVariantEnrichmentAtChromosomalBand**

**1_RegionVsChrBandCount.R**

Count of ASE variants across chromosomal bands and brain regions.

**2_SampleVarSignificance.R**

Permutating ASE variants within samples of each brain region to calculate the P-value assessing whether a chromosomal band was enriched by ASE variants among AD samples or AsymADs or Controls.

**3_HighFreqOrDiffChrBandPlot.R**

Plot for (1) chromosomal bands with highest ASE variant fraction, (2) chromosomal bands exhibiting different patterns of ASE variants enrichment between AD samples and Controls.

<br /><br /><br />




## **ASEVariantAssociationWithClinical**

**1_Boxplot.R**

Analysis and plot for Figure 3A, where comparision of ASE variant frequency between samples grouped by their difference in terms of risk factors (age of death, APOE genotype, and sex) for AD samples and Controls across brain regions.

**2_Heatmap.R**

Analysis and plot for Figure 3B, where fraction of ASE variants showing suggestive or significant associations in the comparison analysis from Figure 3A.

<br /><br /><br />




## **ADAssociatedASEVariantIdentification**

Command of bcftools to extract information of variant of interests from gnomad dataset.

**1_LinearMixedModel.R**

Variant count of predisposing variants in the matched gnomAD ancestry (European of gnomAD is the union of FIN and NFE populations). TCGA population-specific NC P/LPs, exclusively found in a specific TCGA ancestry, are shown as a triangle. Top NC P/LP or top TCGA ancestry-specific NC P/LP, ranked by allele counts in TCGA or gnomAD, was labelled.

**2_ADAssociatedASEUsingOutputFromLinearMFixedModel.R**

(Significance of) Correlations of variant frequencies in the matched ancestries between TCGA and gnomAD.
