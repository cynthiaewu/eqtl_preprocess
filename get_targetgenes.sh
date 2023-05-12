#!/bin/bash

set -e

input_geno=$1
input_results=$2 #result matrix-eQTL file before n_sorted
xqtl_results=$3
top_n=$4
output=$5

var=`sort -rgk 4,4 "${xqtl_results}" | grep -v 'nan' | head -n "${top_n}" | awk -F'\t' '{print $1}'`
for i in $var;
do
    for fn in `find "${input_geno}" -type f -name "genotype_*"`;
    do
        if grep -q -n "$i" "$fn"; then
            filenum=`basename "$fn" | cut -f2 -d'_'`
            head -n 1 "${input_results}_cisregress_${filenum}_sorted" > "${output}/${i}_targetgenes";
            grep $i "${input_results}_cisregress_${filenum}_sorted" | awk -F"\t" '$5<0.05' | sort -gk 5,5 >> "${output}/${i}_targetgenes";
        fi
    done
    # $filenum=`find "${input_geno}" -type f -name "genotype_*" | xargs grep -n $i | cut -f1 -d':' | cut -f2 -d'_'`;
    # echo $filenum
    # head -n 1 "${output}_cisregress_${filenum}_sorted" > "${output}/${var}_targetgenes";
    # grep $var "${output}_cisregress_${filenum}_sorted" | awk -F"\t" '$5<0.05' | sort -gk 5,5 >> "${output}/${var}_targetgenes";
done
