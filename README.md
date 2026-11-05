# Genome Wide and Rare Variants Association Studies Practical

## Introduction:

A Genome Wide Association Study (GWAS) aims to identify genetic markers, often Single Nucleotide Polymorphisms (SNPs), that are associated with a phenotypic trait of interest, for example a disease. The rationale behind a GWAS depends on the type of trait under study:

 - For *binary* phenotypes: the allele frequencies of each SNP are compared between individuals harbouring the trait those who do not. Since we are assuming a random distribution of alleles in the population, large differences in frequencies may indicate an association between the SNP and the trait.

 - For *continuous* phenotypes: the distribution of trait values is compared among individuals carrying different alleles at each SNP. Assuming alleles are randomly distributed in the population, systematic differences in the average trait value between genotypes may indicate an association between the SNP and the trait.

In both cases, a **p-value** is obtained for each SNP, representing the strength of the association between the SNP and the trait. A lower p-value indicates a stronger association. Please note that a low p-value does not imply a large **effect size**, which measures the impact of the SNP on the trait. Effect size is represented as an *odds ratio* for binary traits and as a *beta value* for continuous traits. Note that the effect size also indicates the direction of the allele's impact on the trait (eg: protective of deleterious).

This practical contains three parts:
 - [Quality Control on the input data (QC)](1_QC/)
 - [Genome Wide Association Studies (GWAS)](2_GWAS/)
 - [Rare Variants Association Studies (RVAS)](3_RVAS/)

## The Dataset:

In this practical, we will perform GWAS and RVAS analyses using [PLINK](https://www.cog-genomics.org/plink/1.9/) and R, on a [dataset](./data/) consisting of 90 Asian HapMap individuals (Han Chinese from Beijing and Japanese from Tokyo) from the actual HapMap data. The dataset contains a subset of ~250,000 SNPs from the Affymetrix 500K SNP Array product and a simulated disease phenotype (data adapted from: https://www.cog-genomics.org/plink/1.9/resources).

SNPs were added *in-silico* to specific genes to create a signal for the Rare Variants Association Study (RVAS).

## Starting point:

We are starting with the [Quality Control steps](1_QC/).

## Further Resources:

The GWAS Catalog: https://www.ebi.ac.uk/gwas/ a large repository of summary statistics from published GWAS.

A non-comprehensive list of interesting papers regarding GWAS and RVAS analyses:

### Reviews: 
 - Brandes N., Weissbrod O. & Linial M. (2022) Open problems in human trait genetics. *Genome Biol*
 - Visscher PM. Et al. (2017) 10 Years of GWAS Discovery: Biology, Function, and Translation. *Am J Hum Genet.*
 - Ward LD, Kellis M. (2012) Interpreting noncoding genetic variation in complex traits and human disease. *Nat Biotechnol.*

### Beyond just common SNPs:
 - Harris L. et al. (2025) Genome-wide association testing beyond SNPs. *Nat Rev Genet*
 - Margaux L.A., Hujoel et al. (2022) Influences of rare copy-number variation on human complex traits. *Cell*
 - Auer PL, Lettre G. (2015) Rare variant association studies: considerations, challenges and opportunities. *Genome Med.*

### Working with chrX:
 - Bellavance J et al. (2024) Eight quick tips for including chromosome X in genome-wide association studies. *PLoS Comput Biol.*
 - Le Borgne J. et al. (2025) X‐chromosome-wide association study for Alzheimer’s disease. *Mol Psychiatry*

### Polygenic Risk Scores:
 - Torkamani, A. (2018) The personal and clinical utility of polygenic risk scores. *Nat Rev Genet*
 - Martin, A.R. et al. (2019) Clinical use of current polygenic risk scores may exacerbate health disparities. *Nat Genet*
