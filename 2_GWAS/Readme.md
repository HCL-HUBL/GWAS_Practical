# GWAS (Genome Wide Association Study)

Now that we applied quality control on our dataset, we can move on to the association tests. We are going to perform one association test for each variant present inside our dataset. To analyze if a variant has an impact on the probability to harbor the studied phenotypic trait. 

In a GWAS our null hypothesis (H0) is: `The genetic variant is not associated with the trait`.

## A simple GWAS model:

The model we are going to use depends on the phenotype (Y) under study:
 - For continuous Y: we will use a **linear** regression
 - For categorical Y: we will use a **logistic** regression

> [!Note]
> Other models are possible, for example a model that takes into account dominant or recessive gene action.

Here, we have a case/control phenotype, so we are going to perform a logistic regression on our data. Use the following PLINK command:

```shell
plink --bfile gwas_final --logistic --out gwas_logistic
``` 

This should generate a `gwas_logistic.assoc.logistic` file, containing the results of our regression model. The documentation for this file is available here: https://www.cog-genomics.org/plink/1.9/formats#assoc

**Q2.1:** In the "TEST" column, what does "ADD" means ? (hint: look at the documentation about the linear/logistic model https://www.cog-genomics.org/plink/1.9/assoc#linear)

<details>
    <summary> Answer </summary> 

> 'ADD' means that our model considers an 'additive' effect of the SNP. In other words, we assume that the allelic dosage is additive: homozygous variants have twice the effect of heterozygous variants on the phenotype. The variants are coded 0/1/2 for homozygous reference / heterozygous / homozygous alternative.

</details>

Next, we will explore our results using `R`:

```R
gwas_results <- read.table("gwas_logistic.assoc.logistic", header = TRUE)
gwas_results <- gwas_results[order(gwas_results$P),] # Ordering the results by P-values
```

**Q2.2**: Before looking at the results. How reliable do you think our model is? (hint: go back to the QC part and check if anything could bias our results, or lead to false positive hits). Is there anything missing in our formula?

<details>
    <summary> Answer </summary> 

> During the quality control steps, we noticed that we had two distinct groups in the genetic PCA, with an imbalance in their phenotype distribution. This will likely lead to false positive associations in our GWAS: eg. SNPs that are related to the population groups and not directly to our phenotype of interest. To solve this, we can add the components from the PCA as covariates in our model to account for population structure.

</details>

**Q2.3**: Have a look at the first lines of the data.frame, what are the most associated SNPs? To confirm our hypothesis from the previous questions, make a quick search on gnomAD.

<details>
    <summary> Answer </summary> 

> Our top hit, rs2513514 has an AF of 0.48 in Japanese but only 0.31 in Han Chinese. 
> The 2nd hit rs1548299 has an AF of 0.89 in Japanese but only 0.71 in Han Chinese. 
> Again for the 3rd hit: rs2508756 has an AF of 0.50 in Japanese, and 0.33 in Chinese. 
> Could our results be biased by ancestry?

</details>

## The adjusted model:

Luckily, we have a file in our dataset indicating which population each sample belongs to: `./data/gwas.covars`. This is often not the case, and instead we could use the components of the genetic PCA we computed earlier `pca.eigenvec` as covariates in our model to account for population structure.

Build the model including population as a covariate:

```shell
plink --bfile gwas_final --logistic --covar data/gwas.covars --covar-name POP --out gwas_logistic_adjusted
```

> [!WARNING]
> It is very important to include relevant covariates in your GWAS analysis. Not doing so could result in apparently significant but spurious associations.

As before, a `gwas_logistic_adjusted.assoc.logistic` file is generated, containing the output from our model.

**Q2.4**: Look at the first lines of the file, what is the main difference with the previous one?

<details>
    <summary> Answer </summary> 

> We now have lines with a "POP" value for the "TEST" column. This corresponds to the p-value of the covariate against the phenotype.

</details>

Now that we have working model we can check our top results:

```R
adj_gwas_results <- read.table("gwas_logistic_adjusted.assoc.logistic", header = TRUE)
adj_gwas_results <- adj_gwas_results[adj_gwas_results$TEST == "ADD",] # Only keeping the SNP p-values
adj_gwas_results <- adj_gwas_results[order(adj_gwas_results$P),] # Ordering the results by P-values
```

**Q2.5**: Check the top three SNPs from the initial model: rs2513514, rs1548299 and rs2508756. What are their p-values in the adjusted model? 

<details>
    <summary> Answer </summary> 

> rs2513514 7.76e-04 vs 6.29e-06
> rs1548299 1.90e-02 vs 8.75e-06
> rs2508756 1.02e-03 vs 1.67e-05

</details>

**Q2.6**: Check the top SNPs from the new model in gnomAD.

<details>
    <summary> Answer </summary> 

> rs11204005 has an AF of 0.59 in Japanese, 0.58 in Han Chinese
> rs2460338 has an AF of 0.52 in Japanese, 0.60 in Han Chinese
> rs1052586 has an AF of 0.56 in Japanese, 0.57 in Han Chinese
> These SNPs are not likely to be incluenced by population structure and could be true positive hits!

</details>

Now, we can make the Manhattan and Quantile-Quantile plots:

```shell
head -n 1 gwas_logistic_adjusted.assoc.logistic > gwas_adjusted_SNPs.assoc.logistic # We keep the header line
awk '$5 == "ADD" { print }' gwas_logistic_adjusted.assoc.logistic >> gwas_adjusted_SNPs.assoc.logistic # We only keep the lines corresponding to "ADD"
```

## Manhattan Plot:

```shell
Rscript ../scripts/manhattan_plot.R -i gwas_adjusted_SNPs.assoc.logistic -b BP -m SNP -p P -o plots/manhattan_adjusted.pdf -v TRUE
```

## QQplot:

```shell
Rscript ../scripts/QQplot.R -i gwas_adjusted_SNPs.assoc.logistic -p P -o plots/QQplot_adjusted.png -v TRUE
```

## Post-GWAS steps:

Good! We now have regions significantly associated with our trait. However, we can never be certain that our results are not due to some biases in our cohort / study. After a GWAS it is common to:
 - Validate the results using an independent cohort: with less individuals, but testing only the significant regions. Since we are doing less tests, we do not need to correct as much. If our locus is still significant, it is robust evidence of the association between our locus and the trait.
 - Doing a regional plot: plotting a zoomed version of our significant locus. We can check the LD pattern between the SNPs in that region and look at the neighboring genes.
 - Experimentally confirm the effect of the gene on our trait (With CRISPR knockout for example).

## Pitfalls and other caveats:

### The "missing heritability":

Traits, especially complex ones, are influenced by many factors, including genetics but also environmental. GWAS only captures the genetic part of the trait heritability.

Moreover, GWAS, by design, can only capture common SNPs. But a lot of the genetic components of traits are due to rare, but more impactful SNPs, and larger events (Structural Variants, and Copy Number Variations).

### The "additive" model:

The GWAS model assume a simplistic additive effect of variants. But we know that:
 1. Some variants have dominant/recessive associations with the trait
 2. Some variants have pleiotropic effects (they affect multiple traits at once)
 3. There can be an epistatic effect of different genes (G) for the same trait (G*G interaction)
 4. Some SNPs (G) have an effect only when some environmental (E) factors are in play (G*E)

Some of these issues are easily solved (you can code the SNPs as 0/1 to make a dominant SNP model), others less so.

In the [next section](../3_RVAS/), we will perform Rare Variants Association Studies and try to solve some of the "missing heritability" problem.
