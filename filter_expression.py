import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
import networkx as nx
import random
import copy

def find_iqr(x):
      return np.subtract(*np.percentile(x, [75, 25]))

def normalize_quantiles(df):
    """
    Quantile normalization to the average empirical distribution
    Note: replicates behavior of R function normalize.quantiles from library("preprocessCore")
    Reference:
     [1] Bolstad et al., Bioinformatics 19(2), pp. 185-193, 2003
    Adapted from https://github.com/andrewdyates/quantile_normalize
    """
    M = df.values.copy()

    Q = M.argsort(axis=0)
    m,n = M.shape

    # compute quantile vector
    quantiles = np.zeros(m)
    for i in range(n):
        quantiles += M[Q[:,i],i]
    quantiles = quantiles / n

    for i in range(n):
        # Get equivalence classes; unique values == 0
        dupes = np.zeros(m, dtype=np.int)
        for j in range(m-1):
            if M[Q[j,i],i]==M[Q[j+1,i],i]:
                dupes[j+1] = dupes[j]+1

        # Replace column with quantile ranks
        M[Q[:,i],i] = quantiles

        # Average together equivalence classes
        j = m-1
        while j >= 0:
            if dupes[j] == 0:
                j -= 1
            else:
                idxs = Q[j-dupes[j]:j+1,i]
                M[idxs,i] = np.median(M[idxs,i])
                j -= 1 + dupes[j]
        assert j == -1

    return pd.DataFrame(M, index=df.index, columns=df.columns)


def inverse_normal_transform(M):
    """
    Transform rows to a standard normal distribution
    """
    R = stats.mstats.rankdata(M, axis=1)  # ties are averaged
    if isinstance(M, pd.DataFrame):
        Q = pd.DataFrame(stats.norm.ppf(R/(M.shape[1]+1)), index=M.index, columns=M.columns)
    else:
        Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q

def filter_expr(fname, output, prefix):
    expr = pd.read_csv(fname, sep='\t', index_col=3 )
    expr = expr.loc[:, ~expr.columns.isin(['#chr', 'start', 'end'])]
    num_genes = expr.shape[0]
    log = open(f'{prefix}.log', 'a')

    #print(f'Starting with {num_genes} genes')
    #log.write(f'Starting with {num_genes} genes\n')
    

    #calculate mean, variance, dispersion
    mean = []
    variance = []
    median = []
    max_tpm = []
    for index, gene in expr.iterrows():
        mean.append(np.mean(gene))
        variance.append(np.var(gene))
        median.append(np.median(gene))
        max_tpm.append(max(gene))
        
    gene_info = [list(expr.index), mean, variance, median, max_tpm]
    gene_info = pd.DataFrame(gene_info).T
    gene_info.columns = ['gene', 'mean', 'variance', 'median', 'max_tpm']
    gene_info["disp"] = gene_info.apply(lambda x: x["variance"]/x["mean"], 1)
    
    #calculate iqr
    iqr = []
    for gene, values in expr.iterrows():
        iqr.append(find_iqr(values))
    gene_info.insert(4, 'iqr', iqr)

    gene_info_filter = gene_info[(gene_info['variance']>0) & (gene_info['variance']<50000) & (gene_info['iqr']>0) & (gene_info['median']>0) & (gene_info['max_tpm']>2)]
    expr_filter = expr.loc[gene_info_filter['gene']]
    num_filtered_genes = expr_filter.shape[0]
    print(f'Filtered {num_genes-num_filtered_genes} genes with variance > 0, variance < 50000, iqr < 0, median > 0, and max_tpm > 2')
    log.write(f'Filtered {num_genes-num_filtered_genes} genes with variance > 0, variance < 50000, iqr < 0, median > 0, and max_tpm > 2\n')
    num_genes = num_filtered_genes

    #keep only protein coding genes and filter out AABR0 genes
    rat_gtf = pd.read_csv('/storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/rn6_proteincoding_noAABR0.gtf', sep='\t', header=None)
    protein_genes = []
    for row in rat_gtf[8]:
    #     protein_genes.append(row.split('gene_name "')[1].split('";')[0])
        protein_genes.append(row.split('gene_id "')[1].split('";')[0])
    expr_filter = expr_filter[expr_filter.index.isin(protein_genes)]
    rat_gtf.insert(5, 'gene_id', protein_genes)
    num_filtered_genes = expr_filter.shape[0]
    print(f'Filtered {num_genes-num_filtered_genes} non-protein coding genes')
    log.write(f'Filtered {num_genes-num_filtered_genes} non-protein coding genes\n')
    num_genes = num_filtered_genes
   
    '''
    #filter out genes that fall in seg dup regions
    segdups = pd.read_csv('/gymreklab-tscc/cynthiawu/RatGTEx/Rerun011123/rn6_segdups.bed', sep='\t', header=None)
    segdup_genes = []
    rat_gtf_filtered = rat_gtf[rat_gtf.gene_id.isin(expr_filter.index)]

    for chrom in range(1,21):
        rat_gtf_filtered_chrom = rat_gtf_filtered[rat_gtf_filtered[0]==str(chrom)].sort_values(by=[3])
        segdup_chrom = segdups[segdups[0]==f'chr{chrom}'].sort_values(by=[1])

        for index, row in rat_gtf_filtered_chrom.iterrows():
        #     protein_genes.append(row.split('gene_name "')[1].split('";')[0])
        #     print(row)
            start = row[3]
            end = row[4]
            segdup_chrom_index = (segdup_chrom[1].searchsorted(start))
        #     print(segdup_chrom.iloc[segdup_chrom_index])
        #     print(start)
            #for segdup start position that is before start gene position, if segdup end position is before start gene position
            if segdup_chrom_index > 0:
                if segdup_chrom.iloc[segdup_chrom_index-1][2] > start: 
                    segdup_genes.append(row['gene_id'])
            #for segdup start position that is after start gene position, if segdup start position is before end gene position
            if segdup_chrom_index < len(segdup_chrom):
                if segdup_chrom.iloc[segdup_chrom_index][1] < end:
                    segdup_genes.append(row['gene_id'])
        #     print(segdup_chrom.iloc[segdup_chrom_index-1][2], start)
        #     print(segdup_chrom.iloc[segdup_chrom_index][1], end)
        #     break
        #print(chrom, len(segdup_genes))

    expr_filter = expr_filter.drop(segdup_genes)
    num_filtered_genes = expr_filter.shape[0]
    print(f'Filtered {num_genes-num_filtered_genes} genes in seg dup regions')
    log.write(f'Filtered {num_genes-num_filtered_genes} genes in seg dup regions\n')
    num_genes = num_filtered_genes
    '''

    #calculate gene correlation and filter correlated genes
    expr_corr = expr_filter.T.corr(method='spearman')
    np.fill_diagonal(expr_corr.values, np.nan)
    s = expr_corr.abs().unstack()
    so = s.sort_values(kind="quicksort")
    corr_genes = so[so>0.99]
    
    G = nx.Graph()
    for row in corr_genes.items():
        G.add_edge(row[0][0], row[0][1])    
    keep = []
    remove = set()
    for comp in nx.connected_components(G):
        comp = copy.copy(comp)
        chosen = min(comp)
        keep.append(chosen)
        comp.remove(chosen)
        remove.update(comp)
    
    expr_filter = expr_filter.drop(remove)
    num_filtered_genes = expr_filter.shape[0]
    print(f'Filtered {num_genes-num_filtered_genes} highly correlated genes')
    log.write(f'Filtered {num_genes-num_filtered_genes} highly correlated genes\n')
    num_genes = num_filtered_genes

    #filter out genes with high heterozygousity
    hetrates = pd.read_csv('/storage/cynthiawu/trans_eQTL/RatGTEx/Rerun011123/Scripts/rat-gtex-analysis/preprocessing/blacklist_regions/gene_het_filter', sep='\t')
    hetrates["gene.len"] = hetrates["end"]-hetrates["start"]
    hetrates["snps.per.bp"] = hetrates.apply(lambda x: x["numsnps"]*1.0/x["gene.len"], 1)
    hetrates["hets.per.bp"] = hetrates.apply(lambda x: x["numhets"]*1.0/x["gene.len"], 1)
    hetrates_filtered_expr = hetrates[hetrates['gene'].isin(list(expr_filter.index))]
    total_het = 0
    total_len = 0
    for index, row in hetrates_filtered_expr.iterrows():
        total_het += row['numhets']
        total_len += row['gene.len']
    overall_p = total_het/total_len
    all_p = []
    for index, row in hetrates_filtered_expr.iterrows():
        all_p.append(stats.binom_test(row['numhets'], n=row['gene.len'], p=overall_p, alternative='greater'))
    all_p_corrected = [min(i * len(hetrates_filtered_expr), 1) for i in all_p]
    hetrates_filtered_expr.insert(12, 'p-val', all_p)
    hetrates_filtered_expr.insert(13, 'p-val_corrected', all_p_corrected)
    hetrates_filtered_expr_passgenes = hetrates_filtered_expr[hetrates_filtered_expr['p-val_corrected'] >= 0.05]
    expr_filter = expr_filter[expr_filter.index.isin(list(hetrates_filtered_expr_passgenes['gene']))]
    num_filtered_genes = expr_filter.shape[0]
    print(f'Filtered {num_genes-num_filtered_genes} genes with high heterozygosity')
    log.write(f'Filtered {num_genes-num_filtered_genes} genes with high heterozygosity\n')
   
    num_genes = num_filtered_genes
    print(f'Left with {num_genes} genes')
    log.write(f'Left with {num_genes} genes')

    expr_filter_normalize = normalize_quantiles(expr_filter)
    expr_filter_normalize_inverse = inverse_normal_transform(expr_filter_normalize)
    expr_filter_normalize_inverse.to_csv(output, sep='\t')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fname", required=True, help="Input expression file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered, inverse, normalized expression file")
    parser.add_argument("-p", "--prefix", required=True, help="Expression file prefix")
    params = parser.parse_args()
    filter_expr(params.fname, params.output, params.prefix)


if __name__ == "__main__":
    main()
