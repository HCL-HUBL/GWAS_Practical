# RVAS (Rare Variants Association Studies)

In the last section, we found some variants associated with our trait of interest. However, GWAS, by design, can only find "common" SNPs in our population, as it needs enough carriers to reach statistical significance. We know that rare variants can have a high impact on diseases. For example, a **Mendelian trait** is a trait influenced by a single gene. A known example of Mendelian disease is Huntington disease, often caused by an inherited mutation in the *HHT* gene.

> [!NOTE] 
> The definition of a rare variant differs between studies. Most of the time, variants with a Minor Allele Frequency (MAF) < 0.01 are considered rare.

**Q3.1**: look at the Minor Allele Frequencies (MAF) of the top associated SNPs in our GWAS analysis, are there any rare variants? (Here, we will consider a MAF < 0.01).

```shell
plink --bfile gwas_final --freq
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

Here, we will use 

## Preparation of the input data:

We will start from the quality controlled dataset, just before we removed variants with a MAF < 0.01.