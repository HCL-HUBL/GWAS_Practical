# RVAS (Rare Variants Association Studies)

In the last section, we found some variants associated with our trait of interest. However GWAS, by design, can only find "common" SNPs in our population, as it needs enough carriers of each allele to reach statistical significance. We know that rare variants can have a high impact on diseases. For example, a **Mendelian trait** is a trait influenced by a single gene. A known example of Mendelian disease is Huntington disease, often caused by an inherited mutation in the *HHT* gene.

> [!NOTE] 
> The definition of a rare variant differs between studies. Most of the time, variants with a Minor Allele Frequency (MAF) < 0.01 are considered rare.

**Q3.1**: look at the Minor Allele Frequencies (MAF) of the top associated SNPs in our GWAS analysis, are there any rare variants? (Here, we will consider a MAF < 0.01).

```shell
plink --bfile QC/gwas_final --freq --out RVAS/gwas_final
```

<details>
    <summary> Answer </summary> 

> There are none. Remember that we removed rare SNPs during the last step of the quality control pipeline, to perform less tests. But even if we kept them, it is unlikely that they would reach statistical significance.

</details>

How can we study rare variants, which have low frequencies and are more numerous that common variants?

The main rationale is to aggregate the rare variants by regions (often genes), to aggregate the effect the rare variants located in each region. Multiple tools have been developped to do that, with different approaches.

| Approach | Rationale | Advantage | Disadvantage | Tools |
| --- | --- | --- | --- | --- |
| Burden Test | Collapse rare variants in a genetic region | Powerful if all variants in the region are causal and in the same direction | Less powerful if only a few variants are causal, or have an opposite effect | CAST; WST |
| Adaptative Burden Test | Burden with flexible weights/thresholds | Allows for different directions of effect and for variants-genes interactions | Only applicable in specific cases | KBAC; aSum |
| Variance-component test | Test the variance of genetic effects | Allows for SNP-SNP interaction; Powerful when non-causal variants are present (or protective and deleterious variants) | Less powerful than Burden when all variants are causal and in the same direction | Sequence kernel association test (SKAT) |
| Combined tests | Combines Burden and Variance-component tests | Powerful even with few causal variants, having opposite directions | Needs accurate classification of the variants impact; Choosing an appropriate kernel is often challenging | SKAT-O; MiST |

In this practical, we are going to use PLINK to perform a Burden Test on our dataset. We will consider rare variants as variants with a MAF < 0.01

## Preparation of the input data:

We will start from the quality controlled dataset, just before we removed variants with a MAF < 0.01. This time we will remove variants with a MAF > 0.01, [--max-maf](https://www.cog-genomics.org/plink/1.9/filter#maf) imposes a upper MAF bound for variants.

```shell
plink --bfile QC/gwas_filtered --max-maf 0.01 --make-bed --out RVAS/rvas_final
```

**Q3.2**: Use the same command as above to get the allele frequencies of the `rvas_final` dataset. What is the highest fequency ?

<details>
    <summary> Answer </summary> 

> The highest frequency is now 0.008197. We succesfully removed "common" variants from the dataset.

</details>

## Filtering on chromosomes:

For the purpose of this practical, we are going to focus on chromosomes 8 and 11. We can use the `--chr` parameter in PLINK to only keep variants pertaining to one chromosome:

```shell
plink --bfile RVAS/rvas_final --chr 8,11 --make-bed --out RVAS/rvas_final_subset
```

## Creation of the set test

In PLINK, burden tests are performed by grouping variants into predefined sets (such as genes), collapsing them into a single genetic score per individual, and testing whether this score is associated with the phenotype using linear or logistic regression.

Here, we will define each gene as a set. The coordinates of each gene in the hg19 reference genome are availabile in `./data/glist-hg19`.

```shell
plink --bfile RVAS/rvas_final_subset --make-set data/glist-hg19 --write-set --out RVAS/glist-hg19
```

This will create a `RVAS/glist-hg19.set` file, you can check its content with `less`.

## The association Test:

PLINK collapses variants within each set into a single score per individual.

> [!Note]
> This is a simple burden test, it will only detect associations if all variants affect the trait in the same direction. This lose power if variants in a gene have opposite effects, or many variants are neutral. Other models and tools are available for more complex associations, see the Table at the beginning of this section for more details.

```shell
plink --bfile RVAS/rvas_final_subset --logistic mperm=10000 --assoc set-test perm --set-p 1 --set RVAS/glist-hg19.set --out RVAS/rvas_final_burden
```

**Q3.3**: What is the meaning of `set-test` ?

<details>
    <summary> Answer </summary> 

`set-test` tests the significance of variant sets; https://www.cog-genomics.org/plink/1.9/assoc#set

</details>

**Q3.4**: Whe did we add `--logistic` ?

<details>
    <summary> Answer </summary> 

We have a binary trait (case/controls), so we need to perform a logistic regression to test our sets against the phenotype. If we had a continuous phenotype, we would use `--linear`.

</details>

## Interpretation of results

The output file `rvas_final_burden` should contain the following data:

| SET\_ID | N\_VARIANTS | BETA | P\_VALUE | METHOD |
| --- | --- | --- | --- | --- |
| GENE1 | 10 | 0.5 | 0.001 | burden |
| GENE2 | 8 | \-0.3 | 0.04 | SKAT |

As for the gwas results, we have a `p-value` indicating the strength of association between the set and the phenotype and a `BETA` value indicating the set effect size.

**Q3.5**: Can we interpret the p-value directly ? Is there something we need to do before ? (hint: this is something that applies to GWAS as well).

<details>
    <summary> Answer </summary> 

We are doing one test per set. As there are X sets, we need to correct the p-values for multiple testing. This can be done by adjusting the p-values obtained here, or by using the `--adjust` parameter when performing the test (see below).

</details>

## The adjusted association test:

```shell
plink --bfile RVAS/rvas_final_subset --logistic --adjust --assoc set-test perm --set RVAS/glist-hg19.set --out RVAS/rvas_final_burden_adjusted
```

**Q3.6**: Are there sets (genes) significantly associated with our trait of interest ? How many, and do they make biological sense ?

<details>
    <summary> Answer </summary> 

Yes, ...

</details>

### QQPlot:



## Going further

`test per annotation (SAIGE)`

