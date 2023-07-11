library("MatrixEQTL");

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-g", "--genotype", help="Input genotype file", required=TRUE)
parser$add_argument("-e", "--expression", help="Input expression file", required=TRUE)
parser$add_argument("-o", "--output", help="Output file", required=TRUE)
parser$add_argument("-c", "--covariates", help="Input covariates file")
parser$add_argument("-a", "--snploc", help="SNP location file")
parser$add_argument("-b", "--geneloc", help="Gene location file")
args <- parser$parse_args()
SNP_file_name <- args$genotype
expression_file_name <- args$expression
covariates_file_name <- args$covariates
snploc_file_name <- args$snploc
geneloc_file_name <- args$geneloc
output_file_name_cis <- paste(args$output, "_cis", sep = "")
output_file_name_tra <- paste(args$output, "_trans",  sep = "")

useModel = modelLINEAR;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 5e-2;
pvOutputThreshold_tra = 5e-2;
#pvOutputThreshold = 1;

# Distance for local gene-SNP pairs
cisDist = 5e6;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
#snps$fileOmitCharacters = "-1"; # denote missing values;
snps$fileOmitCharacters = "."; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# Quantile normalization
for( sl in 1:length(gene) ) {
    mat = gene[[sl]];
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(gene) + 1));
    gene[[sl]] = mat;
}
rm(sl, mat);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;
if(length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name);
}

snpspos = read.table(snploc_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(geneloc_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = cvrt,
      output_file_name     = output_file_name_tra,
      pvOutputThreshold     = pvOutputThreshold_tra,
      useModel = useModel,
      errorCovariance = errorCovariance,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos,
      genepos = genepos,
      cisDist = cisDist,
      pvalue.hist = TRUE,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = TRUE);
