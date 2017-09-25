
#Input Files
#GTF="/home/vmason/LinuxShare/gtf/merged.gtf"
GFF="/data2/users/pacholewska/vmason/gff/NCBI/ref_EquCab2.0_top_level.chr.gff3"
BAMINDIR="/data2/users/pacholewska/vmason/bam/HDE9/"
SampleInfoDir="/data2/users/pacholewska/vmason/bam/sampleinfo/RAO_database.VCM.v4.HDE9.no5594.csv"

#Output Files
FEATUREFILE="/data2/users/pacholewska/vmason/counts/features.NCBI.HDE9.tx.tab"
COUNTSOUT="/data2/users/pacholewska/vmason/counts/rawcounts.NCBI.HDE9.tx.tab"
COUNTSOUTTRIM="/data2/users/pacholewska/vmason/counts/rawcounts.NCBI.HDE9.tx.trim.tab"
FPKMFLE="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeq.NCBI.HDE9.FPKM.tx.tab"
OUTFLE="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeqNorm.NCBI.HDE9.NoVSD.tx.tab"
OUTFLEVSD="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeqNorm.NCBI.HDE9.wVSD.tx.tab"
OUTFLERLD="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeqNorm.NCBI.HDE9.wRLD.tx.tab"
OUTFLELG2="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeqNorm.NCBI.HDE9.wLog2.tx.tab"
SF="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/SizeFactors.HDE9.tx.tab"
DISP="/data2/users/pacholewska/vmason/counts/DESeq2Normalized/DESeqDispersion.HDE9.tx.tab"
VARSTABLE="/data2/users/pacholewska/vmason/plots/VSD.HDE9.tx.pdf"
RAWCOUNTHIST="/data2/users/pacholewska/vmason/plots/rawcount.hist.HDE9.tx.pdf"
TRIMCOUNTHIST="/data2/users/pacholewska/vmason/plots/trimcount.hist.HDE9.tx.pdf"
NORMCOUNTHIST="/data2/users/pacholewska/vmason/plots/normcount.hist.HDE9.tx.pdf"
SFHIST="/data2/users/pacholewska/vmason/plots/sizefactors.hist.HDE9.tx.pdf"
VSDCOUNTHIST="/data2/users/pacholewska/vmason/plots/vsdcount.hist.HDE9.tx.pdf"
RLDCOUNTHIST="/data2/users/pacholewska/vmason/plots/rldcount.hist.HDE9.tx.pdf"
LOG2="/data2/users/pacholewska/vmason/plots/log.HDE9.tx.pdf"
RLOG="/data2/users/pacholewska/vmason/plots/rlog.HDE9.tx.pdf"
PCAVSD="/data2/users/pacholewska/vmason/plots/PCAwVSD.HDE9.tx.pdf"
PCARAW="/data2/users/pacholewska/vmason/plots/PCAwRawCounts.HDE9.tx.pdf" # this is not transformed and therefore represents absolute differences not relative differences
PCANORM="/data2/users/pacholewska/vmason/plots/PCAwNormalizedCounts.HDE9.tx.pdf" # this is not transformed and therefore represents absolute differences not relative differences
PCALOG2="/data2/users/pacholewska/vmason/plots/PCAwLOG2.HDE9.tx.pdf"
PCARLOG="/data2/users/pacholewska/vmason/plots/PCAwRLOG.HDE9.tx.pdf"
DISPERSION="/data2/users/pacholewska/vmason/plots/dispersion.HDE9.tx.pdf"
EXPRSAMPLEDIST="/data2/users/pacholewska/vmason/plots/ExpressionSampleDistanceHeatmap.HDE9.tx.pdf"


library( "GenomicFeatures" )
txdb <- makeTxDbFromGFF( GFF, format="gff3", dataSource="NCBI .gff3 rename chr", organism="Equus caballus") # gather features
TxByGene <- transcriptsBy( txdb, by="gene" ) # reduce features to only transcriptsBy( *allfeatures*, by="gene" ) or exonsBy( *allfeatures*, by="gene" )

write.table(TxByGene, file=FEATUREFILE, sep="\t") # write features to file

fls <- list.files( BAMINDIR, pattern="bam$", full=TRUE ) # collect filenames for all bam files in BAMINDIR

library( "Rsamtools" )
bamLst <- BamFileList( fls, yieldSize=1000000 ) # process bam files in chunks #increase number to increase speed and RAM usage, decrease to do the opposite
library(BiocParallel)
param = SnowParam(workers = 8) # specify number of threads # BiocParallel::register(MulticoreParam(workers=5))

library( "GenomicAlignments" )
se <- summarizeOverlaps( TxByGene, bamLst, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE, BPPARAM = param ) # count reads mapped to features
write.table(assay(se), file=COUNTSOUT, sep="\t") # write raw counts to file

library( "DESeq2" )

sampleInfo <- read.csv( SampleInfoDir ) # read in sampleInfo file
sampleInfo <- DataFrame( sampleInfo ) # make sampleInfo a dataframe
sampleInfo$RNAseq_condition <- as.factor(sampleInfo$RNAseq_condition)
sampleInfo$RNAseq_genotype <- as.factor(sampleInfo$RNAseq_genotype)
write("Finished reading in sample info", stderr())


#create index (numbers representing elements) of sampleInfo$BamFile that are in order of colnames(se)
seIdx <- match(colnames(se), sampleInfo$BamFile)
#actually combine the data
colData(se) <- cbind( colData(se), sampleInfo[ seIdx, ] )

dds <- DESeqDataSet( se, design = ~ RNAseq_condition) # RNAseq_genotype + Sex +

pdf(RAWCOUNTHIST)
hist(as.vector(assay(dds)))
dev.off()

#dds <- dds[ rowSums(counts(dds)) > 1, ] # light filtration of rows. Keep row if the sum of the number of mapped reads to feature is > 1
dds <- dds[ rowSums(counts(dds)==0)<(0.1*ncol(assay(se))) ] # replace count matrix with filtered matrix # Keep all rows (features) if the feature has less than 10% 0 raw-count values
write.table(assay(dds), file=COUNTSOUTTRIM, sep="\t") # write filtered raw coutns to file

fpkm <- fpkm(dds, robust=TRUE)
write.table(fpkm, file=FPKMFLE, sep="\t")

pdf(TRIMCOUNTHIST)
hist(as.vector(assay(dds)))
dev.off()

dds <- estimateSizeFactors(dds) # calculate size factors to normalize for library size across samples
normalizedcounts <- counts(dds, normalized=TRUE) # divides all raw counts by size factor, returns the normalized expression matrix
write.table(normalizedcounts, file=OUTFLE, sep="\t") # write DESeq Normalized counts to file
write.table(sizeFactors(dds), file=SF, sep="\t") # write DESeq Size factor estimates for each library to file
dds <- estimateDispersions(dds) # calculate Dispersion
write.table(sizeFactors(dds), file=DISP, sep="\t") # write DESeq dispersion estimates for each gene to file

pdf(NORMCOUNTHIST)
hist(normalizedcounts)
dev.off()
pdf(SFHIST)
hist(sizeFactors(dds))
dev.off()

vsd <- varianceStabilizingTransformation(dds, blind=FALSE) # Transform the dataset to remove dependence of variance on mean # could also use rlog(), but do not use log2()
rld <- rlog(dds, blind=FALSE)

pdf(VSDCOUNTHIST)
hist(as.vector(assay(vsd)))
dev.off()
pdf(RLDCOUNTHIST)
hist(as.vector(assay(rld)))
dev.off()

write.table(assay(vsd), file=OUTFLEVSD, sep="\t") # write DESeq Normalized and variance stabalized counts to file
write.table(assay(rld), file=OUTFLERLD, sep="\t") # write DESeq Normalized and rlog transformed counts to file

notAllZero <- (rowSums(counts(dds))>=0) # Plots all rows where the sum of raw counts is >= 0, Right now this plots everything with the previous filtration if this value is set to 0
log2transf <- dds # renaming variable for log2 transformed data
assay(log2transf) <- log2(counts(dds,normalized=TRUE)[notAllZero,] + 1) # replace raw counts with log2 transformed counts
class(log2transf) <- "DESeqTransform" # Let DESeq2 know that the data has been transformed, and therefore will allow it to be plotted with plotPCA

write.table(assay(log2transf), file=OUTFLELG2, sep="\t") # write DESeq Normalized and rlog transformed counts to file

#FPKM <- fpkm(dds, robust=TRUE)

library(vsn)
pdf(file=VARSTABLE)
meanSdPlot(assay(vsd[notAllZero,])) # plot variance stabalized counts
dev.off()
pdf(file=LOG2)
meanSdPlot(assay(log2transf))
dev.off()
pdf(file=RLOG)
meanSdPlot(assay(rld[notAllZero,]))
dev.off()
pdf(file=PCAVSD) # plot PCA of VSD transformed expression count values
plotPCA(vsd, intgroup=c("RNAseq_condition", "RNAseq_genotype")) # colored by family
dev.off()
rawdds <- dds
class(rawdds) <- "DESeqTransform" # trick DESeq2 into thinking these coutns have been transformed
pdf(file=PCARAW) # plot PCA of RLOG transformed expression count values
plotPCA(rawdds, intgroup=c("RNAseq_condition", "RNAseq_genotype")) # colored by family
dev.off()
normdds <- dds
assay(normdds) <- normalizedcounts # replace raw counts with normalized counts
class(normdds) <- "DESeqTransform" # trick DESeq2 into thinking these coutns have been transformed
pdf(file=PCANORM) # plot PCA of normalized count values ( not transformed )
plotPCA(normdds, intgroup=c("RNAseq_condition", "RNAseq_genotype")) # colored by family
dev.off()
pdf(file=PCALOG2) # plot PCA of Log2() transformed expression count values
plotPCA(log2transf, intgroup=c("RNAseq_condition", "RNAseq_genotype")) # colored by family
dev.off()
pdf(file=PCARLOG) # plot PCA of RLOG transformed expression count values
plotPCA(rld, intgroup=c("RNAseq_condition", "RNAseq_genotype")) # colored by family
dev.off()
pdf(file=DISPERSION)
plotDispEsts(dds)
dev.off()

sampleDists <- dist( t( assay(rld) ) )
library("pheatmap")
library("RColorBrewer")
pdf(file=EXPRSAMPLEDIST) # FAILED
sampleDistMatrix <- as.matrix( sampleDists )

rownames(sampleDistMatrix) <- paste( colData(vsd)$RNAseq_condition, rownames(colData(vsd)), sep="-" ) # make names of labels for pheatmap
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
sampleDists <- dist(t(assay(vsd)))
dev.off()

save.image(file="/data2/users/pacholewska/vmason/logs/HDE9.RData")
savehistory(file="/data2/users/pacholewska/vmason/logs/HDE9.Rhistory")