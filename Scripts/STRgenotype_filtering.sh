#!/bin/bash

set -e

input_folder=$1
output_folder=$2

#printf  "Starting STR filtering" > "${output_folder}/STR_filtering.log"
#for i in {1..9} {11..20}; 
#    do 
#    var_count0=`grep -v '^#' "${input_folder}/chr${i}.vcf" | wc -l`
#    printf  "Starting with ${var_count0} variants in chr${i}\n" >> "${output_folder}/STR_filtering.log"
#    bcftools view -f "PASS,."  "${input_folder}/chr${i}.vcf" > "${output_folder}/chr${i}_filter_callrate.vcf";
#    var_count1=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate.vcf" | wc -l`
#    printf  "Filtered $((var_count0-var_count1)) variants that did not pass 0.8 call rate\n" >> "${output_folder}/STR_filtering.log"
#    dumpSTR --min-locus-hwep 10e-6 --min-locus-het 0.1 --use-length --vcf  "${output_folder}/chr${i}_filter_callrate.vcf" --out "${output_folder}/chr${i}_filterformat_callrate_het";
#    bcftools view -f "PASS,." "${output_folder}/chr${i}_filterformat_callrate_het.vcf" > "${output_folder}/chr${i}_filter_callrate_het.vcf";
#    var_count2=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het.vcf" | wc -l`
#    printf  "Filtered $((var_count1-var_count2)) variants that did not pass heterozygousity 0.1 and hwe 10e-6\n" >> "${output_folder}/STR_filtering.log"
#    bedtools intersect -a "${output_folder}/chr${i}_filter_callrate_het.vcf" -b data/rn6_segdups.bed -v -header > "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf";
#    var_count3=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf" | wc -l`
#    printf  "Filtered $((var_count2-var_count3)) variants that are in segmental duplication regions\n" >> "${output_folder}/STR_filtering.log"
#    #bedtools intersect -a  "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf" -b data/rn6_exon_noAABR0_TSS3kb.bed -u -header > "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf"; 
#    bedtools intersect -a  "${output_folder}/chr${i}_filter_callrate_het_nosegdups.vcf" -b data/rn6_keep_snp_regions.bed -u -header > "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf"; 
#    var_count4=`grep -v '^#' "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf" | wc -l`
#    printf  "Filtered $((var_count3-var_count4)) variants that are not in exon or TSS+/-3kb regions\n" >> "${output_folder}/STR_filtering.log"
#    vcftools --vcf "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon.vcf" --extract-FORMAT-info GB --out "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon_genotype";
#    printf  "Left with  $((var_count4)) variants in chr${i}\n" >> "${output_folder}/STR_filtering.log"
#    echo "Finished chr${i}";
#done
#
#cp "${output_folder}/chr1_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT"
#for i in {2..20}; do tail -n +2 "${output_folder}/chr${i}_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" >> "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" || continue; done
#
#python formatSTRgenotype.py -g "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype.GB.FORMAT" -o "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype"

python get_unique_geno.py -i "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype" -o "${output_folder}/allchr_filter_callrate_het_nosegdups_proteincoding_exon_genotype_unique"
