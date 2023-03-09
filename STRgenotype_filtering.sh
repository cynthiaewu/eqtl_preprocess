#!/bin/bash

set -e

input_folder=$1
output_folder=$2
#geno_label=$2

#var_count0=`cat $input_vcf | grep -v "^#" | wc -l | cut -f1 -d' '`
#printf  "Starting with ${var_count0} variants\n" > "${geno_label}.log"

printf  "Starting STR filtering" > "${output_folder}/STR_filtering.log"
for i in {1..9} {11..20}; 
    do 
    var_count0=`grep -v '^#' "${input_folder}/chr${i}.vcf" | wc -l`
    printf  "Starting with ${var_count0} variants in chr${i}\n" >> "${output_folder}/STR_filtering.log"
    bcftools view -f "PASS,."  "${input_folder}/chr${i}.vcf" > "${output_folder}/chr${i}_filter_callrate.vcf";
    var_count1=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate.vcf" | wc -l`
    printf  "Filtered $((var_count0-var_count1)) variants that did not pass 0.8 call rate\n" >> "${output_folder}/STR_filtering.log"
    dumpSTR --min-locus-het 0.1 --use-length --vcf  "${output_folder}/chr${i}_filter_callrate.vcf" --out "${output_folder}/chr${i}_filterformat_callrate_het"; 
    bcftools view -f "PASS,." "${output_folder}/chr${i}_filterformat_callrate_het.vcf" > "${output_folder}/chr${i}_filter_callrate_het.vcf";
    var_count2=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het.vcf" | wc -l`
    printf  "Filtered $((var_count1-var_count2)) variants that did not pass heterozygousity 0.1\n" >> "${output_folder}/STR_filtering.log"
    bedtools intersect -a "${output_folder}/chr${i}_filter_callrate_het.vcf" -b data/rn6_segdups.bed -v -header > "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf";
    var_count3=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf" | wc -l`
    printf  "Filtered $((var_count2-var_count3)) variants that are in segmental duplication regions\n" >> "${output_folder}/STR_filtering.log"
    bedtools intersect -a  "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf" -b data/rn6_exon_noAABR0_TSS3kb.bed -u -header > "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf"; 
    var_count4=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf" | wc -l`
    printf  "Filtered $((var_count3-var_count4)) variants that are not in exon or TSS+/-3kb regions\n" >> "${output_folder}/STR_filtering.log"
    vcftools --vcf "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf" --extract-FORMAT-info GB --out "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon_genotype";
    printf  "Left with  $((var_count4)) variants in chr${i}\n" >> "${output_folder}/STR_filtering.log"
    echo "Finished chr${i}";
done

cp "${output_folder}/chr1_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT"
for i in {2..20}; do tail -n +2 "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" >> "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" || continue; done

python formatSTRgenotype.py -g "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" -o "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype"


#  for i in {1..20}; do dumpSTR --min-locus-het 0.1 --vcf STR_genotypes_rn6/chr$i\_filter_callrate.vcf --out STR_genotypes_rn6/chr$i\_filterformat_callrate_het; done
#  for i in {1..20}; do bcftools view -f "PASS,."  STR_genotypes_rn6/chr$i\_filterformat_callrate_het.vcf > STR_genotypes_rn6/chr$i\_filter_callrate_het.vcf; done
##  for i in {1..20}; do bedtools intersect -a STR_genotypes_rn6/chr$i\_filter_callrate_het.vcf -b /storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/Scripts/eqtl_preprocess/data/rn6_segdups.bed -v -header > STR_genotypes_rn6/chr$i\_filter_callrate_het_nosegdups.vcf; done
#   for i in {1..20}; do bedtools intersect -a STR_genotypes_rn6/chr$i\_filter_callrate_het_nosegdups.vcf -b /storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/Scripts/eqtl_preprocess/data/rn6_exon_noAABR0_TSS3kb.bed -u -header > STR_genotypes_rn6/chr$i\_filter_callrate_het_nosegdups_proteincoding_exon.vcf; done
#
#   for i in {2..20}; do vcftools --vcf STR_genotypes_rn6/chr$i\_filter_callrate_het_nosegdups_proteincoding_exon.vcf --extract-FORMAT-info GB --out STR_genotypes_rn6/chr$i\_filter_callrate_het_nosegdups_proteincoding_exon_genotype; done
#   for i in {2..20}; do tail -n +2 chr2_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT >> allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT; done
#
# for i in {1..20}; do dumpSTR --min-locus-het 0.1 --vcf /gymreklab-tscc/eveloff/rn6_calling_pipeline/dumpstr_out/chr$i\.vcf --out STR_genotypes_rn6/chr$i\_filter_het; done
#
#'''
##segdups
##bcftools annotate --rename-chrs data/chrom_map --threads 8 -o "${geno_label}_filter.recode.rename.vcf.gz" -Oz "${geno_label}_filter.recode.vcf"
#
#bedtools intersect -a "${geno_label}_filter.recode.rename.vcf.gz" -b data/rn6_segdups.bed -v -header > "${geno_label}_filter.recode.rename.nosegdups.vcf"
#
#var_count1=`wc -l "${geno_label}_filter.recode.vcf" | cut -f1 -d' '`
#
#var_count2=`wc -l "${geno_label}_filter.recode.rename.nosegdups.vcf" | cut -f1 -d' '`
#printf  "Filtered $((var_count1-var_count2)) variants in segdup regions\n" >> "${geno_label}.log"
#
##TSS +/-3kb and exons of protein coding genes only (-u to report at most once for each entry)
#bedtools intersect -a "${geno_label}_filter.recode.rename.nosegdups.vcf" -b data/rn6_exon_noAABR0_TSS3kb.bed -u -header > "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf"
#
#var_count3=`wc -l "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" | cut -f1 -d' '`
#printf  "Filtered $((var_count2-var_count3)) variants not in TSS +/-3kb or exons of protein coding genes\n" >> "${geno_label}.log"
#
##LD pruning
##plink --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" --indep-pairwise 200 100 0.99 --recode vcf --out "${geno_label}_filter_nosegdups_proteincoding_exon"
#
##plink --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.vcf" --extract "${geno_label}_filter_nosegdups_proteincoding_exon.prune.in"  --recode vcf --out "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned"
#
##var_count4=`wc -l "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" | cut -f1 -d' '`
#printf  "LD pruned out $((var_count3-var_count4)) variants\n" >> "${geno_label}.log"
#
#var_count5=`grep -v "^#" "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" | wc -l | cut -f1 -d' '`
#printf  "Left with ${var_count5} variants\n" >> "${geno_label}.log"
#
##convert into 012 genotype format
#vcftools --vcf "${geno_label}_filter.recode.rename.nosegdups.proteincoding_exon.pruned.vcf" --012 --out "${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr"
#
#python convert012_matrix.py -f "${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr.012" -o "${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr_genotype"
#
#sed -i -E '1 s/([A-Z0-9]*)_//g' "${geno_label}_filter_nosegdups_proteincoding_exon_pruned_allchr_genotype"
#''' 
