SampleInfoDir="D:\\LinuxShare\\Programs\\DESeq2\\DataIN\\HDE9\\RAO_database.VCM.v4.HDE9.no559443.csv"
RAWREADCOUNTS="D:\\LinuxShare\\Programs\\DESeq2\\DataIN\\HDE9\\rawcounts.NCBI.HDE9.tx.sort.tab"
FEATUREFLE="D:\\LinuxShare\\Programs\\DESeq2\\DataIN\\HDE\\features.NCBI.alldose.tx.trim.tab"

FPKMFLE="\\home\\vmason\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\DESeq.NCBI.FPKM.tx.tab"
NORMFLE="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\DESeqNorm.NCBI.NoVSD.RSum1Trim.HDE9.tx.tab"
NORMFLEm24="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\DESeqNorm.NCBI.NoVSD.NoTrim.HDE9.m24.tx.tab"
SF="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\SizeFactors.HDE9.m24.tx.tab"
NORMCOUNTHIST="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\Plots\\normcount.hist.HDE9.m24.tx.pdf"
SFHIST="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\Plots\\sizefactors.hist.HDE9.m24.tx.pdf"
VSDFLE="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\DESeq.NCBI.HDE9.m24.wVSD.tx.tab"
VSDHIST="D:\\LinuxShare\\Programs\\DESeq2\\NormalizedCounts\\Plots\\VSD.hist.HDE9.m24.tx.pdf"

library( "DESeq2" )
sampleInfo <- read.csv( SampleInfoDir ) # read in sampleInfo file
sampleInfo <- DataFrame( sampleInfo ) # make sampleInfo a dataframe
sampleInfo$Mare <- as.factor(sampleInfo$Mare)
sampleInfo$Age <- as.numeric(sampleInfo$Age)
sampleInfo$Fam1 <- as.factor(sampleInfo$Fam1)
sampleInfo$Fam2 <- as.factor(sampleInfo$Fam2)
sampleInfo$RNAseq_condition <- as.factor(sampleInfo$RNAseq_condition) # don't need to because binary quantitative variables are equivalent to binary categorical (dummy) variables 

countdata <- read.table(RAWREADCOUNTS, sep="\t", header=TRUE)

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleInfo, design = ~ 1 + Mare + Age + Fam1 + Fam2 + RNAseq_condition ) # By default, the functions in DESeq2 package will use the last variable in the formula for building result tables and plotting

dim(counts(dds))
dds <- dds[ rowSums(counts(dds)) > 1, ] # require each gene to have at least one read ailgned
#dds <- dds[ rowSums(counts(dds)==0)<(0.1*ncol(assay(se))) ] # require 90% of individuals to have >= 1 count
dim(counts(dds))
#dds <- dds[ rowMeans(counts(dds)) >= 24, ] # subset(counts(dds), rowMeans(counts(dds)) >= 24) # require row (gene) mean read count to be above or equal to 24. # 24 was determined through KS Test on untrimmed normalized RNAseq coutns.
#dim(counts(dds))
dds <- estimateSizeFactors(dds) # calculate size factors to normalize for library size across samples
normalizedcounts <- counts(dds, normalized=TRUE) # divides all raw counts by size factor, returns the normalized expression matrix
write.table(normalizedcounts, file=NORMFLE, sep="\t", quote=FALSE) # write DESeq Normalized counts to file
write.table(sizeFactors(dds), file=SF, sep="\t", quote=FALSE) # write DESeq Size factor estimates for each library to file
dds <- dds[ rowMeans(counts(dds, normalized=TRUE)) >= 24, ]
dim(counts(dds))

normalizedcounts <- counts(dds, normalized=TRUE) # divides all raw counts by size factor, returns the normalized expression matrix
write.table(normalizedcounts, file=NORMFLE, sep="\t", quote=FALSE) # write DESeq Normalized counts to file

res <- results(dds, contrast=c(RNAseq_condition))

pdf(NORMCOUNTHIST)
hist(normalizedcounts)
dev.off()
pdf(SFHIST)
hist(sizeFactors(dds))
dev.off()

dds <- dds[ rowSums(counts(dds)==0)<(0.1*ncol(countdata)) ] 

dds <- varianceStabilizingTransformation(dds, blind=FALSE) # blind = FALSE to utilize the groups defined in experimental design

write.table(assay(dds), file=VSDFLE, sep="\t", quote=FALSE) # write DESeq Normalized counts to file
pdf(VSDHIST)
hist(as.vector(assay(dds)))
dev.off()


