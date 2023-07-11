import argparse
import pandas as pd

def add_genotypes(x):
    if x == '.':
        return '.'
    geno = x.split('|')
    return int(geno[0]) + int(geno[1])

def intersectFiles(geno_file, geno_out):
    STR_data = pd.read_csv(geno_file, sep='\t')
    STR_data.index = STR_data['CHROM'] + '_' + STR_data['POS'].astype(str)
    STR_data = STR_data.drop(['CHROM', 'POS'], axis=1)
    STR_data = STR_data.applymap(add_genotypes)
    #STR_data = STR_data[~STR_data.eq(STR_data.iloc[:, 0], axis=0).all(1)]
    num_var = len(STR_data)
    print(f'Left with {num_var} STRs')
    STR_data.to_csv(geno_out, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotype", required=True, help="Input STR genotype file")
    parser.add_argument("-o", "--geno_out", required=True, help="Output genotype file of intersecting samples")
    params = parser.parse_args()
    intersectFiles(params.genotype, params.geno_out)


if __name__ == "__main__":
    main()
