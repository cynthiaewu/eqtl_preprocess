#!/bin/bash

set -e

genotype=$1
expression=$2
covariates=$3
output=$4
snploc=$5
geneloc=$6
genotype_segments=$7 #folder to split the genotype file into

echo "Running Matrix-eqtl to get cis-eqtls"
#Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm_cistrans.r -g $genotype -e $expression -c $covariates -o $output -a $snploc -b $geneloc

echo "Regressing out cis eqtls from expression"
#python cis_regress.py -e $expression -g $genotype -c "${output}_cis" -o "${expression}_cisregress"

#split -l 500 -d $genotype "${genotype_segments}/genotype_"

echo "Running Matrix-eqtl on expression with cis-eqtls regress out"

for file in "${genotype_segments}/*";
do
    n=$(echo $file | sed -E 's/^.*_//');
    #echo "$file"; echo "$file $file";
    #echo "${output}_cisregress${n}"
    Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_cisregress" -c $covariates -o "${output}_cisregress_${n}"
done

