import argparse
import pandas as pd
import sklearn.decomposition

def perform_pca(input, num_factors, output):
    #expression = pd.read_csv(input, sep='\t', index_col=0)
    expression = pd.read_csv(input, index_col=0)
    samples = list(expression.columns)
    expression_trans = expression.values.transpose()
    pca = sklearn.decomposition.PCA(n_components=num_factors)
    pca.fit(expression_trans)
    expression_trans_pca = pca.transform(expression_trans)
    PCs = ['PC' + str(i) for i in range(len(expression_trans_pca[0]))]
    PC_trans_df = pd.DataFrame(expression_trans_pca.T, columns=samples, index=PCs)
    PC_trans_df.to_csv(output, sep='\t')
    #return pca, expression_trans_pca, PC_trans_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input expression file")
    parser.add_argument("-n", "--num_factors", required=True, type=int, help="Number of PCs")
    parser.add_argument("-o", "--output", required=True, help="Output PCs file")
    params = parser.parse_args()
    perform_pca(params.input, params.num_factors, params.output)


if __name__ == "__main__":
    main()
