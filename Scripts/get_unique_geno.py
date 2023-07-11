import argparse
import pandas as pd
from collections import Counter

def get_uniquegeno(fname, out):
    df = pd.read_csv(fname, sep='\t', index_col=0)
    exclude = []
    for index, row in df.iterrows():
        counts = (Counter(row))
        nunique = len(counts)

        # at least 3 unique genotypes (not counting '.' missing)
        if '.' in counts:
            nunique = nunique - 1
        if nunique >= 3:
            ngeno = 0
            raregeno = []
            for element in counts:
                if element != '.':
                    #at least 3 samples per genotype
                    if counts[element] >= 3:
                        ngeno += 1
                    else:
                        raregeno.append(element)
            # if at least 3 samples for at least 3 genotypes, check for rare genotypes (less than 3 samples) and replace those as missing
            if ngeno >= 3:
                if raregeno:
                    df.loc[index] = df.loc[index].replace(raregeno, ['.']*len(raregeno))
            else:
                exclude.append(index)
        else:
            exclude.append(index)
    df = df.drop(index=exclude)
    df.to_csv(out, sep='\t')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fname", required=True, help="Input genotype file")
    parser.add_argument("-o", "--out", required=True, help="Output genotype file with unique genotypes with at least 3 samples each")
    params = parser.parse_args()
    get_uniquegeno(params.fname, params.out)


if __name__ == "__main__":
    main()
