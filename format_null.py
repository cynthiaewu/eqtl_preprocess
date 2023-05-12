import argparse
import pandas as pd

def get_nullfile(output_name, null_out, count):
    all_snps_run = pd.read_csv(f'{output_name}_00_sorted_xqtl/results.tab', sep='\t')
    for i in range(1,9):
        cur_df = pd.read_csv(f'{output_name}_0{i}_sorted_xqtl/results.tab', sep='\t')
        all_snps_run = pd.concat([all_snps_run, cur_df], axis=0)
    for i in range(10,count):
        cur_df = pd.read_csv(f'{output_name}_{i}_sorted_xqtl/results.tab', sep='\t')
        all_snps_run = pd.concat([all_snps_run, cur_df], axis=0)
    all_snps_run = all_snps_run.dropna()
    all_snps_run[['CPMA', 'xQTL']].to_csv(null_out, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_name", required=True, help="Input output file name")
    parser.add_argument("-n", "--null_out", required=True, help="Output empirical null file with xQTL and CPMA values")
    parser.add_argument("-c", "--count", required=True, type=int, help="Number of files that should be merged")
    params = parser.parse_args()
    get_nullfile(params.output_name, params.null_out, params.count)


if __name__ == "__main__":
    main()
