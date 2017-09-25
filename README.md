# HorseeQTLs2017

Code was ran with R v.3.3, Python 2.7.13, and GNU bash, version 4.1.2(2)-release (x86_64-redhat-linux-gnu)
Code is provided in directory ./HorseeQTLs2017/scripts/
Code was executed in the order specified in the flow chart provided in ./HorseeQTLs2017/images/FigS1Supplement.2.1.tif

#1 Rscript summarizeOverlapsToDESeqTwo.6.0.tx.r
1 Requires R libraries library( "GenomicFeatures" ), library( "Rsamtools" ), library(BiocParallel), library( "GenomicAlignments" ), library( "DESeq2" ), library(vsn), library("pheatmap"), library("RColorBrewer")
1.a Input: path to directory with indexed .bam files
1.b Input: path to .gff3 annotation file for the aligned genome
1.c Input: path to sample information file (comma delimited)
1.d By default attempts to utilize 8 threads. To change this modify line 48 ==> param = SnowParam(workers = 8)

#2 python KSTest.VCM.py
2 Requires Python 2.7, numpy, pandas, matplotlib, scipy.stats, pylab
2.a This code is modified to accept a gene feature count matrix from source code by dmnfarrell: https://github.com/dmnfarrell/mirnaseq
2.a Input: Gene expression count matrix, individuals as columns rows are genes, and tab delimited.
2.a Input: First header line should only have individual names, and no geneid header.
2.a id1 id2 id3
2.a MC1R 5 20 14
2.a ASIP 10 13 8

#3 Rscript NormalizeAndTrimRawCounts.tx.2.R
3 Requires DESeq2
3.a Requires path to feature file (including covariates)
3.b Requires path to read count file (tab delimited)
2.b Input: Gene expression count matrix, individuals as columns rows are genes, and tab delimited.
2.b Input: First header line should only have individual names, and no geneid header.
3.c Requires path to sample information file (comma delimited)
3.a Modify line 37 to represent the KS test cutoff ==> dds <- dds[ rowMeans(counts(dds, normalized=TRUE)) >= 24, ]

#4 python PEER.OutputFactorsAndResiduals.wCov.HDE9.py
4 Requires PEER
4 Requires python modules matplotlib, matplotlib.backends.backend_pdf, pylab, os
4.a Input: current working directory
4.b Input: KS test trimmed, normalized, variance stabalized gene expression count matrix
4.b Format: Genes are rows and Individuals as columns, and tab delimited. First header line should only have individual names, and no geneid header.
4.c Input: covariate file
4.c Format: Matrix is individuals as rows and covariates as columns. Header line should have one name per column

#Call Tag SNPs
#5 python BatchFastTaggerPerChromsome.py
5 NOTE: must have separated input files by chromosome (with samtools: tabix SNPs.vcf chr1 > SNPs.chr1.vcf )
5 Requires: python modules os, glob, subprocess, datetime

#6 python ConvertFastTaggerTagSNPsToMatrixeQTLAndeQTLBMA.py
6 Function: gather FastTagger output and reformats to Matrix eQTL and eQTLBMA input files
6 Requires: python modules os, glob, pandas, datetime

#Matrix eQTL
#7 Rscript MatrixeQTL.Linear.r
7 Requires library MatrixEQTL
7 ### Requires individuals in all files to be in the same order in all input files ###
7.a Input: SNP genotype matrix
7.a Format: coded as 0,1,2 (0 = homozygous reference, 1 = heterzygous, 2 = homozygous alternative)
7.a Format: SNPs are rows and Individuals as columns, and tab delimited. ADD 'snpid\t' TO THE BEGINNING OF THE HEADER LINE this is followed by individual names
7.b Input: SNP location file
7.c Input: Gene epression matrix (PEER residuals)
7.c Format: Genes are rows and Individuals as columns, and tab delimited. ADD 'geneid\t' TO THE BEGINNING OF THE HEADER LINE this is followed by individual names
7.d Input: Gene location file
7.e Input (optional): covariate file ## this shouldn't be needed as known and unknown covariates would have been regressed out from gene expression matrix in PEER

#8 python MakeDFsForeQTLStripchartsAndLinearModel.wFDRcut.py
8 NOTE: It might be beneficial to limit the number of eQTLs considered for the following two scripts (8&9). i.e. Limit to only significant eQTLs, or the most significant eQTLs.
8 Function: Makes one tab delimited file for each gene / SNP pair in current directory. Each file will have ordered information for the following script to make linear model plots
8.a Input: path to eQTL output file from Matrix eQTL
8.b Input: path to covariate file input into Matrix eQTL
8.c Input: path to gene epression file input into Matrix eQTL
8.c Input: path to SNP genotype file input into Matrix eQTL

#9 Rscript MakeLinearModelPlots.Linear.wEQTLResiduals.9.3.r
9 Function: Determine high confidence eQTLs by calculating cooks distance and leverage, and test if global assumption is passed with gvlma. Plot eQTL.
9 Requires: R libraries car, ggplot2, plyr, gvlma
9.a Input: path to directory where python MakeDFsForeQTLStripchartsAndLinearModel.wFDRcut.py was ran.

#eQTLBMA
The code below (10 & 11) utilizes files formatted for eQTLBMA
#10 bash eQTLBMA4.1.bf.hm.bfs.EstPie0.sh
10 Function: basically to calculate Pi0

#11 bash eQTLBMA4.2.bfs.bestSNP1.sh
11 Input: input Pi0 value calculated from script 10 before running this code



