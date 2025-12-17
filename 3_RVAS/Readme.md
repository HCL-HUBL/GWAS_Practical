# RVAS (Rare Variants Association Studies)

In the last section, we performed GWAS on our dataset, but only captured "common" SNPs, indeed rare SNPs, by definition, are not present in enough samples to reach statistical significance. 

**Q3.1**: look at the Minor Allele Frequencies (MAF) of the top associated SNPs in our GWAS analysis, are there any rare variants ? 

<details>
    <summary> Answer </summary> 

> There are none. Remember that we even removed rare SNPs during the last step of quality control pipeline, to perform less tests ! But even if some were present, they would not reach statistical significance in the logistic regression model that we used.

</details>

There are 