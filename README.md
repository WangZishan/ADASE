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

Calculate the P-value assessing whether one variant exhibites allele specific expression based on binomial test.

**2_ASEVariantIdentification.R**

ASE detection from bi-allelic heterozygous variants of autosomal chromosomes that pass stringent quality control in each sample.

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

Analysis and plot for Figure 3A, where comparision of ASE variant frequency between samples grouped by their difference in terms of risk clinical factors (age of death, APOE genotype, and sex) for AD samples and Controls across brain regions.

**2_Heatmap.R**

Analysis and plot for Figure 3B, where fraction of ASE variants showing suggestive or significant associations in the comparison analysis from Figure 3A.

<br /><br /><br />




## **ADAssociatedASEVariantIdentification**

**1_LinearMixedModel.R**

Application of linear mixed effects model onto allele expression profile for calculating (BH-corrected) P-value assessing whether one variant exhibites differential allelic imbalance between AD samples and Controls.

**2_ADAssociatedASEUsingOutputFromLinearMixedModel.R**

Identification of AD associated ASE variant using output from the linear mixed effects model as input.
