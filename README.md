# eqtl_preprocess

## Filter SNP genotype with 
* include sites with minor allele frequency>0.1
* include only sites in exons and TSS +/- 3kb of protein coding genes
* exclude sites that has more than 0.8 missing data 
* exclude sites that fall in segmental duplication regions

```
bash genotype_filtering.sh $genotype $geno_label
```
$genotype is the genotype file in 012 format for all SNPs

$geno_label is the label subsequent files are named. e.g. 'IL_LHb_NAcc_OFC_PL'

## Format STR genotypes
1. Format STR genotypes as dosages (repeat numbers added up for both alleles)
```
formatSTRgenotype.py --genotype --geno_out
```
--genotype is the .GB.FORMAT file from plink for STR

--geno_out is the formatted STR genotype file

## Filter STR genotype with 
* include sites with heterozygousity>0.1
* include only sites in exons and TSS +/- 3kb of protein coding genes
* exclude sites that has more than 0.8 missing data 
* exclude sites that fall in segmental duplication regions

```
bash genotype_filtering.sh $input_folder $output_folder
```
$input_folder is the folder that contains all formatted STR genotype files separated by chromosomes

$output_folder is the folder to write output to

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
$expression is the expression file 

$expr_label is the label subsequent files are named. e.g. 'OFC'

## Regressing out cis effects and Matrix-eQTL

1. Run Matrix-eQTL and obtain cis-eQTLs (within 1mb of gene and pval<0.05)
2. Run elastic net model to regress out effects of cis-eQTLs
3. Run Matrix-eQTL on expression dataset with regressed out cis effects
```
bash cisregress_matrixeqtl.sh $genotype $expression $covariates $output $snploc $geneloc $genotype_segments
```
$genotype is the genotype file

$expression is the expression file

$covariates is the covariates file

$output is the location to write output files

$snploc is the file with the locations of variants in format 'snp chr pos'

$geneloc is the file with the locations of genes in format 'geneid chr pos1 pos2'

$genotype_segments is the path to the folder where genotype files are broken into files of 500 variants
