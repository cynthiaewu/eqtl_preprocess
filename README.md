# eqtl_preprocess

## Filter SNP genotype with 
* include sites with minor allele frequency>0.1
* include only sites in exons and TSS +/- 3kb of protein coding genes
* exclude sites that has more than 0.8 missing data 
* exclude sites that fall in segmental duplication regions

```
bash genotype_filtering.sh $genotype $geno_label
```

## Filter STR genotype with 
* include sites with heterozygousity>0.1
* include only sites in exons and TSS +/- 3kb of protein coding genes
* exclude sites that has more than 0.8 missing data 
* exclude sites that fall in segmental duplication regions

```
bash genotype_filtering.sh $genotype $geno_label
```

## Filter and preprocess expression dataset
1. Apply filtering steps
  * include protein coding genes only
  * include genes that have ‘variance'>0, 'iqr'>0
  * exclude genes that fall in segmental duplication regions
  * exclude highly correlated genes. corr_genes>0.99

2. Quantile normalize 
3. Run PEER
```
bash expression_pipeline.sh $expression $expr_label
```

## Regressing out cis effects and Matrix-eQTL

1. Run Matrix-eQTL and obtain cis-eQTLs (within 1mb of gene and pval<0.05)
2. Run elastic net model to regress out effects of cis-eQTLs
3. Run Matrix-eQTL on expression dataset with regressed out cis effects
```
bash cisregress_matrixeqtl.sh $genotype $expression $covariates $output $snploc $geneloc $genotype_segments
```
