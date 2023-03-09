#!/bin/bash

set -e

genotype=$1
expression=$2
covariates=$3
output=$4
snploc=$5
geneloc=$6
genotype_segments=$7 #folder to split the genotype file into

python intersect_samples.py -g $genotype -e $expression -c $covariates -o "${genotype}_intersect" -p "${expression}_intersect" -q "${covariates}_intersect"

echo "Running Matrix-eqtl to get cis-eqtls"
Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm_cistrans.r -g "${genotype}_intersect" -e "${expression}_intersect" -c "${covariates}_intersect" -o $output -a $snploc -b $geneloc

echo "Regressing out cis eqtls from expression"
python cis_regress.py -e "${expression}_intersect" -g "${genotype}_intersect" -c "${output}_cis" -o "${expression}_intersect_cisregress"

genotype_intersect="${genotype}_intersect";
#echo $genotype_intersect
export genotype_intersect
tail -n +2 $genotype_intersect | split -d -l 500 - --filter='sh -c "{ head -n 1 $genotype_intersect; cat; }" > $FILE' "${genotype_segments}/genotype_"

 echo "Running Matrix-eqtl on expression with cis-eqtls regress out"
 for file in $(find ${genotype_segments} -type f);
 do
     n=$(echo $file | sed -E 's/^.*_//');
     Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_intersect_cisregress" -c "${covariates}_intersect" -o "${output}_cisregress_${n}"
     (head -n 1 "${output}_cisregress_${n}" && tail -n +2 "${output}_cisregress_${n}" | sort -gk 1) > "${output}_cisregress_${n}_sorted"
     xQTL-run --input "${output}_cisregress_${n}_sorted" --out "${output}_cisregress_${n}_sorted_xqtl" --cpma --xqtl --null_method chi2 --threads 14  
 
 done
