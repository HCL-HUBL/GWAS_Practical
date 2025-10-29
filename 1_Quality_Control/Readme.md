# Quality Control

In this pratical the quality control steps will be shared between the common single variants (GWAS) and rare variants (RVAS) association studies.

## The data:

The data for this practical is available under `./data/` and contains:
 - *gwas.{bim,bed,fam}*: genotyping data in the [plink](https://www.cog-genomics.org/plink/1.9/) format
 - *gwas.covars*: a file of covariates
 - *regions.bed*: a file listing regions to analyse

### The plink format:

To perform a GWAS we need the phenotype of each individual and their genotype at each marker. To perform the quality control steps, we will use the [PLINK](https://www.cog-genomics.org/plink/1.9/) software, which read the required information from coupled .bim / .bed and .fam files.

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

 1. Have a look in chr15.bim. At which genomic position is located SNP rs4932079 ? What are two the observed alleles for this SNP is this dataset ?

 <details>
    <summary> Answer </summary> 
    rs7111951 is located on chr11 at position 19973135 and the two observed alleles are *G* and *A*.
 </details>

 2. Have a look in gwas.fam. Which gender is individual NA18550 ?

 <details>
    <summary> Answer </summary> 
    NA18547 is a female (5th column = 2)
 </details>

## Our first steps with plink:

PLINK is a command line standalone software dedicated to GWAS analyses. You will call this software with the alias `plink``. Try to print the help menu with the following command-line:

```shell
plink --help
```

For each of its function, PLINK will print a summary of data and actions in the console. Just load the gwas dataset with PLINK by using the `--bfile`` option, indicating you are using a binary dataset :

```shell
plink --bfile gwas --noweb
```

3. How many individuals contains this dataset ? How many males/females ? How many variants ?

 <details>
    <summary> Answer </summary> 
    The dataset contains 90 people (45 males, 45 females) and 261,184 variants. 
 </details>

## Filtering low quality data:


