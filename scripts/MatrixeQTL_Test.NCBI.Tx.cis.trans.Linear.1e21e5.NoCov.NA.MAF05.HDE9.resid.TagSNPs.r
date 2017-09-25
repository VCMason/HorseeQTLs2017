# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL);

## Location of the package with the data files.
base.dir = "/data2/users/pacholewska/vmason/matrixeqtl/";
# base.dir = '.';

#Give information for filenames
model <- "Linear"
normalization <- "DESeqVSD"
treatment <- "HDE9"
pvals <- "c1e2t1e5.TagSNPs"
covs <- "NA"
covar <- "NoErrCov" # no error covariance matrix
sv <- "PEERResid" # noSV == no surragate variables # wSV is with surrogate varaibles ( from SVA ) # wFac is with factors from PEER
na <- "NoNA" # there is no missing data in SNP file
MAF <- "MAF0.05"
feature <- "tx"


## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "dataIN/HDE9_no559443/82_RAO_HD.2M.HDE9.trim2eQTLIndivs.ne1e3.MateQTL.TagSNP.tab", sep="");
snps_location_file_name = paste(base.dir, "dataIN/HDE9_no559443/82_RAO_HD.2M.sort.trim2eQTLIndivs.ne1e3.MateQTL.AllSNPLoc.tab", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "dataIN/HDE9_no559443/residuals.HDE9.VSD.MatrixeQTL.sort.tx.tab", sep="");
gene_location_file_name = paste(base.dir, "dataIN/HDE9_no559443/features.NCBI.HDE9.trim.MatEQTL.tx.tab", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() # paste(base.dir, "dataIN/HDE9_no559443/RAO_database.VCM.v4.HDE9.no5594.Mateqtl.binary.wFactors.csv", sep="");

# Output file name
output_file_name_cis = paste(base.dir, paste("Cis", model, normalization, treatment, pvals, covs, covar, sv, na, MAF, feature, "out", sep="."), sep="");
output_file_name_tra = paste(base.dir, paste("Trans", model, normalization, treatment, pvals, covs, covar, sv, na, MAF, feature, "out", sep="."), sep="");
output_file_name_cis2 = tempfile();
output_file_name_tra2 = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-2; #1e-2
pvOutputThreshold_tra = 1e-5; #1e-5

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
#errorCovariance = read.table(paste(base.dir, "/dataIN/Covariance.AllDose.SNPv6.clean.MAF0.05.NoHead.tab", sep=""), sep="\t");
#write(dim(errorCovariance), stderr());
#errorCovariance <- as.matrix(errorCovariance); # as.matrix coerces the nxn list into nxn numeric matrix (of type double) # as.numeric() would convert list to double (but of only one dimension)
#write(dim(errorCovariance), stderr());
#write(typeof(errorCovariance), stderr());

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header=TRUE, stringsAsFactors = FALSE, sep="\t");
genepos = read.table(gene_location_file_name, header=TRUE, stringsAsFactors = FALSE, sep="\t");

write("Files loaded into memory", stderr())
write(class(genepos), stderr())
write(class(genepos[1,3]), stderr())
write(class(genepos[2,3]), stderr())
write(gene_location_file_name, stderr())
#write(genepos[1,], stderr())
#write(genepos[1,3], stderr())
write(class(genepos[1, 3]) %in% c("integer", "numeric"), stderr())
write(class(genepos[2, 3]) %in% c("integer", "numeric"), stderr())

meh = Matrix_eQTL_main(
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
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

write("First eQTL analysis for histogram finished", stderr())
## Plot the histogram of all p-values
pdf(file = paste(base.dir, paste("histogram.cis.trans", model, normalization, treatment, pvals, covs, covar, sv, na, MAF, feature, "out", "pdf", sep="."), sep=""));

plot(meh);

dev.off();
write("Histogram plotted successfully", stderr())
##

meq = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra2,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis2,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra2);
unlink(output_file_name_cis2);

write("Second eQTL analysis for QQ-Plot finished", stderr())

## Results:

cat('Analysis done in: ', meh$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(meh$cis$eqtls);
cat('Detected distant eQTLs:', '\n');
show(meh$trans$eqtls);

## Plot the QQ-Plot of all p-values

pdf(file = paste(base.dir, paste("QQPlot.cis.trans", model, normalization, treatment, pvals, covs, covar, sv, na, MAF, feature, "out", "pdf", sep="."), sep=""));

plot(meq);

dev.off();
