#!/bin/bash

set -e

input=$1
expr_label=$2

echo "Starting filtering"
zcat $input | sed -e '2,$s/^/chr/' > "${expr_label}.expr.tpm.rename.bed"
gene_count0=`wc -l "${expr_label}.expr.tpm.rename.bed" | cut -f1 -d' '`
printf "Starting with $((gene_count0-1)) genes\n" > "${expr_label}.log"
echo "Starting with $((gene_count0-1)) genes" 

bedtools intersect -a "${expr_label}.expr.tpm.rename.bed" -b /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/rn6_segdups.bed -v -header > "${expr_label}.expr.tpm.rename.nosegdups.bed"  

gene_count1=`wc -l "${expr_label}.expr.tpm.rename.nosegdups.bed" | cut -f1 -d' '`
printf "Filtered $((gene_count0-gene_count1)) genes in segdup regions\n" >> "${expr_label}.log"
echo "Filtered $((gene_count0-gene_count1)) genes in segdup regions"

python /storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/eqtl_preprocess/filter_expression.py -f "${expr_label}.expr.tpm.rename.nosegdups.bed" -o "${expr_label}.expr_filter_normalize_inverse_proteincoding_nosegdups.tsv" -p $expr_label

echo "Finished filtering"
echo "Starting PEER"
path=`pwd`
Rscript /storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/eqtl_preprocess/run_PEER.r "$path" "${expr_label}.expr_filter_normalize_inverse_proteincoding_nosegdups.tsv" "${expr_label}_expr_filter_normalize_inverse_proteincoding_nosegdups"

sed -i -e '1 s/"X//g' -e 's/"//g' "${expr_label}_expr_filter_normalize_inverse_proteincoding_nosegdups_peer_residuals.tsv"
 
