import argparse
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNetCV
import warnings
from sklearn.exceptions import ConvergenceWarning
from collections import defaultdict

def elasticnet(gene_cisqtls_ldprune, expression, genotype):
    notchosen = 0
    index = 0
    num_warnings = 0
    for gene in gene_cisqtls_ldprune:
        if gene in gene in expression.index:
            Y = expression.loc[gene]
            X = []
            for snp in gene_cisqtls_ldprune[gene]:
                X.append(list(genotype.loc[snp]))
            X = np.array(X).T
            regr = ElasticNetCV(cv=8, max_iter=5000)

            try:
                warnings.filterwarnings("error", category=ConvergenceWarning, module="sklearn")
                fit = regr.fit(X, Y)
            except:
                num_warnings += 1
            else:
                X = X.astype(np.float64)
                pred = fit.predict(X)
                residuals = Y - pred
                expression.loc[gene] = residuals

        else: 
            notchosen += 1
    return expression, notchosen, num_warnings


def apply_elasticnet(f_expression, f_genotype, f_cisqtls, output):
    expression = pd.read_csv(f_expression, sep='\t', index_col=0)
    genotype = pd.read_csv(f_genotype, sep='\t', index_col=0)
    cis_qtls = pd.read_csv(f_cisqtls, sep='\t')
    gene_cisqtls = defaultdict(list)
    for index, row in cis_qtls.iterrows():
        gene = row['gene']
        if len(gene_cisqtls[gene]) < 10:
            gene_cisqtls[gene].append(row['SNP'])
    
    expression_cisregress, notchosen, num_warnings = elasticnet(gene_cisqtls, expression, genotype)
    expression_cisregress.to_csv(output, sep='\t')
        


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--f_expression", required=True, help="Input expression file")
    parser.add_argument("-g", "--f_genotype", required=True, help="Input genotype file")
    parser.add_argument("-c", "--f_cisqtls", required=True, help="Input matrixeqtl file with cis-eqtls only")
    parser.add_argument("-o", "--output", required=True, help="Output expression file with cis-eqtls regressed out")
    params = parser.parse_args()
    apply_elasticnet(params.f_expression, params.f_genotype, params.f_cisqtls, params.output)


if __name__ == "__main__":
    main()
