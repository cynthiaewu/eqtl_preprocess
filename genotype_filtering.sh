#!/bin/bash

set -e

input_vcf=$1
geno_label=$2

var_count0=`zcat $input_vcf | grep -v "^#" | wc -l | cut -f1 -d' '`
printf  "Starting with ${var_count0} variants\n" > "${geno_label}.log"

vcftools --gzvcf $input_vcf --out "${geno_label}_filter" --recode --maf 0.1 --max-missing 0.8

#segdups
bcftools annotate --rename-chrs /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/chrom_map --threads 8 -o "${geno_label}_filter.recode.rename.vcf.gz" -Oz "${geno_label}_filter.recode.vcf"

bedtools intersect -a "${geno_label}_filter.recode.rename.vcf.gz" -b /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/rn6_segdups.bed -v -header > "${geno_label}_filter.recode.rename.nosegdups.vcf"

var_count1=`wc -l "${geno_label}_filter.recode.vcf" | cut -f1 -d' '`

var_count2=`wc -l "${geno_label}_filter.recode.rename.nosegdups.vcf" | cut -f1 -d' '`
printf  "Filtered $((var_count1-var_count2)) variants in segdup regions\n" >> "${geno_label}.log"

#TSS +/-3kb and exons of protein coding genes only
bedtools intersect -a "${geno_label}_filter.recode.rename.nosegdups.vcf" -b /storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/rn6_exon_noAABR0_TSS3kb.bed -header > "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf"

var_count3=`wc -l "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" | cut -f1 -d' '`
printf  "Filtered $((var_count2-var_count3)) variants not in TSS +/-3kb or exons of protein coding genes\n" >> "${geno_label}.log"

#LD pruning
plink --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" --indep-pairwise 200 100 0.99 --recode vcf --out "${geno_label}_filter_nosegdups_proteincoding_exon"

plink --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" --extract "${geno_label}_filter_nosegdups_proteincoding_exon.prune.in"  --recode vcf --out "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned"

var_count4=`wc -l "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" | cut -f1 -d' '`
printf  "LD pruned out $((var_count3-var_count4)) variants\n" >> "${geno_label}.log"

var_count5=`grep -v "^#" "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" | wc -l | cut -f1 -d' '`
printf  "Left with ${var_count5} variants\n" >> "${geno_label}.log"

#convert into 012 genotype format
vcftools --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" --012 --out "genotype/${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr"

python /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/convert012_matrix.py -f "genotype/${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr.012" -o "genotype/${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr_genotype"

sed -i -E '1 s/([A-Z0-9]*)_//g' "genotype/${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr_genotype"
 
