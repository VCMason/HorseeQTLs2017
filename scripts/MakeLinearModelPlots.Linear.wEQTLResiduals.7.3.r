library(car)
library(ggplot2)
library(plyr)
library(gvlma)

sigcov = 2 # the position of the covariate in the linear model (remember to add one since R adds the intercept by default if not specified)
cookscut = 0.5 # cooks distance threshold to determine outliers
#INDIR = "/home/vmason/LinuxShare/Programs/MatrixEQTL/Plots"
INDIR = "D:\\LinuxShare\\Programs\\MatrixEQTL\\dataOUT\\MCK1_MAF0.05_VSD_no559443\\Linear_TagSNPswExpressionResiduals\\Cis"

outputfile = paste(INDIR, "Out_CooksDistanceAndLeverage.txt", sep = "\\")
nonoutliersfile = paste(INDIR, "Out_NoOutliers.CooksDistanceAndLeverage.txt", sep = "\\")
nofalsefile = paste(INDIR, "Out_NoOutliers.NoFalseAssumptions.txt", sep = "\\")
globalfile = paste(INDIR, "Out_GlobalStatPass.txt", sep = "\\")
globalNoOutlierfile = paste(INDIR, "Out_GlobalStatPass.NoOutliers.txt", sep = "\\")
globalfailfile = paste(INDIR, "Out_GlobalStatFail.txt", sep = "\\")

fls <- list.files( INDIR, pattern="tab$", full=TRUE ) # collect filenames for all tab files in INDIR

count = 0
count2 = 0
count3 = c(0,0,0,0,0)
count4 = 0
count5 = 0
trueeqtls = c()
globaleqtls = c() # eQTLs that pass the global stat in gvlma
globalNoOutliereqtls = c() # eQTLs that pass the global stat in gvlma and do not have outliers
globalfaileqtls = c() # eQTLs that fail the global stat in gvlma

output = paste('File', 'Cooks.Distance.Outliers', 'Leverage.Outliers', 'Cooks.Distance.Outlier.Individuals', 'Leverage.Outlier.Individuals', sep='\t')
write(output, file=outputfile)
write(output, file=nonoutliersfile)

for (f in fls)
{
	print(f)
	d <- read.table(file=f, sep="\t", header=TRUE)
	rownames(d) <- d[,1]
	
	# multiple linear regression model
	lm.model <- lm(d[,3] ~ d[,2], data=d) # d[,3] == gene expression, d[,2] == snp genotype, d[,4:(ncol(d)-1)] are other covariates besides RNAseq_condition, d[,ncol(d)] is RNAseq_condition
	#lm(d[,3] ~  d[,2] + d[,4] + d[,5] + d[,6] + d[,7] + d[,8] + d[,9] + d[,10] + d[,11] + d[,ncol(d)] + d[,2]:d[,ncol(d)]), data=d)
	#lm(d[,3] ~ d[,2] + d[,4:(ncol(d)-1)] + d[,2]:d[,ncol(d)])
	#lm(d[,3] ~ . - BamFile - Gene, data=d)
	
	cooks <- cooks.distance(lm.model)
	cooks.na <- cooks[is.na(cooks)] # keeps only NA values
	cooks <- cooks[!is.na(cooks)] # keeps only non-NA (i.e. numerical) values # by removing NA values the cookouts filter works properly
	cookouts <- cooks[cooks>cookscut] # inforces cooks distance cutoff to find outliers
	cookouts <- c(cooks.na, cookouts) # joins the NA values as outliers because these are commonly the only individual of a genotype for with a case/control status and can drive significance by themselves.
	
	gate1 = 0
	leverage <- hatvalues(lm.model)
	highmoments <- leverage[leverage > 2*(length(lm.model$coefficients)/length(lm.model$residuals))]
	
	if (length(cookouts)>0 & length(highmoments)>0) # if there is an outlier value for cooks distance, or for leverage. Then report it.
		{
			output = paste(f, paste(cookouts, collapse=','), paste(highmoments, collapse=','), paste(names(cookouts), collapse=','), paste(names(highmoments), collapse=','), sep="\t")
			count = count + 1
			write(output, file=outputfile, append=TRUE)
		}
	else
		{
			gate1 = 1
			out = paste(f, paste(cookouts, collapse=','), paste(highmoments, collapse=','), paste(names(cookouts), collapse=','), paste(names(highmoments), collapse=','), sep="\t")
			count2 = count2 + 1
			write(out, file=nonoutliersfile, append=TRUE)
		}
	assumptions <- gvlma(lm.model)
	capture.output( summary(assumptions), file=paste(f,'gvlma', sep='.') )
	#decisions <- summary(assumptions)$Decision
	decisions <- display.gvlmatests(assumptions)$Decision
	gate2 = 1
	for ( i in 1:5 )
		{
		if ( "Assumptions NOT satisfied!" %in% decisions[i] )
			{
			count3[i] = count3[i] + 1
			gate2 = 0
			} 
		}
	if ( gate2 == 0 ) { count4 = count4 + 1 } else { count5 = count5 + 1 } # count4 if eQTL failed at least one assumption # count5 if eQTL failed not assumptions.
	if ( gate1 == 1 & gate2 == 1 ) { trueeqtls <- c(trueeqtls, f) } # are there not outlier individuals and no false assumptions for this eQTL? if not then append filename.
	if ( decisions[1] == "Assumptions acceptable." ) { globaleqtls <- c(globaleqtls, f) } else { globalfaileqtls <- c(globalfaileqtls, f) }
	if ( decisions[1] == "Assumptions acceptable." & gate1 == 1 ) { globalNoOutliereqtls <- c(globalNoOutliereqtls, f) }
	print(paste(gate1,gate2))
	
	pdf(file=paste(f, "pdf", sep="."))
	
	#dSummary <- d %>% group_by(d[,2]) %>% summarize(m = mean(d[,3]), se = sqrt(var(d[,3])/length(d[,3]))) # requires dplyr # d[,3] == gene expression, d[,2] == snp genotype
	
	simple <- lm(d[,3] ~ d[,2]) # simple linear regression model with genotype as independent variable lm(gene_expression ~ gentoype)
	# d[,3] == gene expression, d[,2] == snp genotype
	p <- ggplot(data=d, aes(x = d[,2], y = d[,3], color = factor(d[,2]))) +
	  geom_violin(trim=TRUE, show.legend=FALSE) + theme_bw() +
	  scale_color_manual(values=c(1, 2, 3)) +
	  geom_point(shape = factor(d[,ncol(d)]), color = as.factor(d[,ncol(d)]+1), position = position_jitter(width = 0.10, height = 0.0), alpha = 0.9) + # need to add 1 to RNAseq_condition because 0 is not a color value
	  geom_smooth(method = "lm", se = TRUE, fill = "gray90", color = "black", formula = y ~ x ) +
	  theme(text=element_text(size=8)) +
	  xlab("Genotype") +
	  ylab("Gene Expression") +
	  ggtitle(paste("Slope = ",signif(simple$coef[[sigcov]], 5))) # paste("Adj R2 =",signif(summary(simple)$adj.r.squared, 5), " Intercept =",signif(simple$coef[[1]],5 ), " P =",signif(summary(simple)$coef[sigcov,4], 5))
	print(p)
	ggsave(paste(f, "Gen.pdf", sep="."), plot=p, width=2.5, height=2.5)
	ggsave(paste(f, "Gen.tiff", sep="."), plot=p, width=2.5, height=2.5, dpi=600)
	# Simple linear regression plot of lm(gene expression ~ gentoype)
	# this block added individual labels # d[,3] == gene expression, d[,2] == snp genotype
	p <- p + geom_text(aes(label = rownames(d)), nudge_x = 0.3, size = 1.75, color="black")
	print(p)
	
	dev.off()
}
print(paste(count, '= number of eQTLs with high cooks disatnce and high leverage'))
print(paste(count2, '= number of eQTLs that do not have both high cooks distance AND high leverage'))
print(count3)
print('= number of eQTLs that fail each assumption (Global Stat, Skewness, Kurtosis, Link Function, Heteroscedasticity) tested by gvlma')
print(paste(count4, '= number of eQTLs that fail at least one assumption tested by gvlma'))
print(paste(count5, '= number of eQTLs that do not fail any assumptions tested by gvlma'))
print(paste(round((count4/(count4+count5))*100, 2), '= percentage of eQTLs that fail at least one assumption tested by gvlma'))
write(paste(trueeqtls, collapse='\n'), nofalsefile)
write(paste(globaleqtls, collapse='\n'), globalfile)
print(paste(length(globalNoOutliereqtls), '= number of eQTLs that pass global stat in gvlma and no outliers with high Cooks Distance and high Leverage'))
write(paste(globalNoOutliereqtls, collapse='\n'), globalNoOutlierfile) # This is the file I would consider high confidence eQTLs (even though 'trueeqtls' is a more stringent filter). # eQTLs pass global stat in gvlma and has no outliers
write(paste(globalfaileqtls, collapse='\n'), globalfailfile)