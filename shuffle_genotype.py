import argparse
import pandas as pd
import random

def intersectFiles(geno_file, geno_out):
    genotype = pd.read_csv(geno_file, sep='\t', index_col=0)
    intersect_order = list(genotype.columns)
    sample_labels = list(genotype.columns)
    random.shuffle(sample_labels)
    genotype.columns = sample_labels
    genotype_shuffled = genotype[intersect_order]
    genotype_shuffled.to_csv(f'{geno_file}_shuffled', sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotype", required=True, help="Input genotype file")
    parser.add_argument("-o", "--geno_out", required=True, help="Output genotype file of shuffled samples")
    params = parser.parse_args()
    intersectFiles(params.genotype, params.geno_out)


if __name__ == "__main__":
    main()
