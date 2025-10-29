# Genome Wide and Rare Variants Association Studies

## Introduction

A Genome Wide Association Study (GWAS) aims at discovering genomic alleles that are associated to a phenotypic trait of interest.
The main rationale behind a GWAS is to compare SNPs frequencies between groups of people harboring the trait or not.

Here, we will perform a GWAS analysis on 90 Asian HapMap individuals (Han Chinese from Beijing and Japanese from Tokyo). From the actual HapMap SNP data, ~250,000 SNPs have been extracted, which are the autosomal SNPs on half of the Affymetrix 500K SNP Array product. Along with a simulated disease phenotype.

This pratical will cover the standard steps of an association analysis: 
 - Quality Control on the input data
 - Single Variants association testing (GWAS)
 - Rare Variants association testing (RVAS)
 - Visualisation and interpretation of the results

## Data:

The data for this practical is available under `./data/` and contains:
 - *gwas.{bim,bed,fam}*: genotyping data in the [plink](https://www.cog-genomics.org/plink/1.9/) format
 - *gwas.covars*: a file of covariates
 - *regions.bed*: a file listing regions to analyse

To perform a GWAS we need the phenotype of each individual and their genotype at each marker. To perform the quality control steps, we will use
the [PLINK](https://www.cog-genomics.org/plink/1.9/) software, which read the required information from coupled .bim / .bed and .fam files.

 - The [.bim](https://www.cog-genomics.org/plink/1.9/formats#bim) file lists all the genotyped variants (position and alleles):
  
    | chromosome | identifier | position\_morgans | position\_bp | allele\_1 | allele\_2 |
    | --- | --- | --- | --- | --- | --- |
    | 1 | rs7519837 | 0 | 1596068 | C | T |
    | 1 | rs3737628 | 0 | 1755094 | T | C |
    | 1 | rs7511905 | 0 | 1825948 | A | C |
    | 1 | rs3855951 | 0 | 1836464 | C | T |
    | 1 | rs6603803 | 0 | 1844850 | G | A |

> [!NOTE]  
> The actual file does **not** have a header file.

 - The [.fam](https://www.cog-genomics.org/plink/1.9/formats#fam) file contains sample information:

| FID | IID | paternal\_ID | maternal\_ID | sex | phenotype |
| --- | --- | --- | --- | --- | --- |
| CH18524 | NA18524 | 0 | 0 | 1 | 1 |
| CH18526 | NA18526 | 0 | 0 | 2 | 1 |
| CH18529 | NA18529 | 0 | 0 | 2 | 1 |
| CH18532 | NA18532 | 0 | 0 | 2 | 1 |
| CH18537 | NA18537 | 0 | 0 | 2 | 2 |

 - The [.bed](https://www.cog-genomics.org/plink/1.9/formats#bed) file represents the genotype calls. Must be accompanied by a .bim and .fam file.

> [!NOTE]  
> More information about the formats can be found in the Plink documentation: https://www.cog-genomics.org/plink/1.9/formats.

**Questions**: 

 - Have a look in chr15.bim. At which genomic position is located SNP rs4932079 ? What are two the observed alleles for this SNP is this dataset ?

 <details>
    <summary> Answer </summary> 
     rs7111951 is located on chr11 at position 19973135 and the two observed alleles are *G* and *A*.
 </details>

## Quality Control

PLINK is a command line standalone software dedicated to GWAS analyses. You will call this software with the alias `plink``. Try to print the help menu with the following command-line:

```shell
plink --help
```
