#!/bin/bash

set -e

genotype=$1
expression=$2
covariates=$3
output=$4
genotype_segments=$5 #folder to split the genotype file into

python intersect_samples.py -g $genotype -e $expression -c $covariates -o "${genotype}_intersect" -p "${expression}_intersect" -q "${covariates}_intersect"

#get shuffled genotype labels
python shuffle_genotype.py -g "${genotype}_intersect" -o "${genotype}_intersect_shuffled"  
#get empirical null
genotype_intersect="${genotype}_intersect_shuffled";
echo $genotype_intersect
export genotype_intersect
mkdir -p ${genotype_segments}_shuffled
tail -n +2 $genotype_intersect | split -d -l 500 - --filter='sh -c "{ head -n 1 $genotype_intersect; cat; }" > $FILE' "${genotype_segments}_shuffled/genotype_"

 echo "Running Matrix-eqtl on expression with cis-eqtls regress out"
 for file in $(find ${genotype_segments}_shuffled -type f);
 do
     n=$(echo $file | sed -E 's/^.*_//');
     #Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_intersect" -c "${covariates}_intersect" -o "${output}_shuffled_${n}"
     Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_intersect" -o "${output}_shuffled_${n}"
     (head -n 1 "${output}_shuffled_${n}" && tail -n +2 "${output}_shuffled_${n}" | sort -gk 1) > "${output}_shuffled_${n}_sorted"
     xQTL-run --input "${output}_shuffled_${n}_sorted" --out "${output}_shuffled_${n}_sorted_xqtl" --cpma --xqtl --null_method chi2 --threads 14  

 done

count=`ls  ${genotype_segments}_shuffled | wc -l`
python format_null.py -o "${output}_shuffled" -n "${output}_shuffled_nullresults" -c $count

#run pipeline on observed genotype data

genotype_intersect="${genotype}_intersect";
#echo $genotype_intersect
export genotype_intersect
mkdir -p ${genotype_segments}
tail -n +2 $genotype_intersect | split -d -l 500 - --filter='sh -c "{ head -n 1 $genotype_intersect; cat; }" > $FILE' "${genotype_segments}/genotype_"

 echo "Running Matrix-eqtl on expression with cis-eqtls regress out"
 for file in $(find ${genotype_segments} -type f);
 do
     n=$(echo $file | sed -E 's/^.*_//');
     #Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_intersect" -c "${covariates}_intersect" -o "${output}_${n}"
     Rscript Run_Matrix_eQTL_PC_PEER_quantile_norm.r -g $file -e "${expression}_intersect" -o "${output}_${n}"
     (head -n 1 "${output}_${n}" && tail -n +2 "${output}_${n}" | sort -gk 1) > "${output}_${n}_sorted"
#     xQTL-run --input "${output}_cisregress_${n}_sorted" --out "${output}_cisregress_${n}_sorted_xqtl" --cpma --xqtl --null_method chi2 --threads 14  
#     xQTL-run --input "${output}_cisregress_${n}_sorted" --out "${output}_cisregress_${n}_sorted_xqtl_nullrun1" --cpma --xqtl --null_method eigen --precomputed_null /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/OFC/shuffled_results_run1/OFC_nullresults_shuffled_run1 --threads 14  
     #echo "${output}_cisregress_${n}_sorted"
     xQTL-run --input "${output}_${n}_sorted" --out "${output}_${n}_sorted_xqtl_nullrun" --cpma --xqtl --null_method eigen --precomputed_null "${output}_shuffled_nullresults" --threads 14  
     #xQTL-run --input "${output}_cisregress_${n}_sorted" --out "${output}_cisregress_${n}_sorted_xqtl_nullrun" --cpma --xqtl --null_method eigen --precomputed_null /gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/OFC/shuffled_results_run1/OFC_nullresults_shuffled_run1  --threads 14  
 done

