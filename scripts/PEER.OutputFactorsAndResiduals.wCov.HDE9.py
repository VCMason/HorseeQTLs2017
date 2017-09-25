''' module add SequenceAnalysis/Filtering/peer/1.3 '''
#python -c "import Module_Name" 
#expr = SP.loadtxt('/home/vmason/MyPrograms/peer-master/examples/data/expression.csv', delimiter=',')
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import pylab as PL
import os
#cwd = '/data2/users/pacholewska/vmason/PEER' # os.getcwd()
cwd = '/data2/users/pacholewska/vmason/PEER/dataIN/HDE9'
infile = "DESeq.NCBI.HDE9.m24.wVSD.tx.2.tab"
covfile = "RAO_database.VCM.v4.HDE9.no559443.CovOnly.csv"

outfile = "residuals.HDE9.VSD.tx.tab"
outfactorfile = 'factors.HDE9.VSD.tx.tab'
pdffile = 'PEER.plots.HDE9.VSD.tx.pdf'

def Plot_Alpha(Alpha, color="blue"):
	fig = PL.figure()
	PL.plot(1.0 / Alpha,lw=4, color=color)
	min_a,max_a = (1.0/Alpha).min(), (1.0/Alpha).max()
	PL.ylim(min_a - 0.1*(max_a - min_a), max_a + 0.1*(max_a - min_a))
	PL.xlabel("Factors")
	PL.ylabel("Factor relevance")
	return(fig)
	
def PlotGraph(data):
	fig = PL.figure()
	PL.plot(data)
	return(fig)
	
def WriteOUT(residuals, header, genenames):
	OUT = open("%s/%s" % (cwd, outfile), "w")
	print residuals.shape
	print residuals
	print type(residuals)
	newres = [genenames] # make the first line the header line of gene names
	for r, g in map(None, residuals, header.strip().split("\t")):
		print len(r)
		print len(header)
		temp = [g] # add the individual ID to the beginning of every row
		for i in r:
			temp.append(str(i)) # make all the float numbers strings so we can join() them
		newres.append(temp)
		
	OUT.write('\n'.join(['\t'.join(r) for r in newres])) # put the puzzle together
	OUT.close()
	return()
	
def WriteOUTFactors(factors, header):
	OUT = open("%s/%s" % (cwd, outfactorfile), "w")
	OUT.write('\n'.join(['\t'.join(map(str, f)) for f in factors])) # same as [[str(f) for f in indiv] for indiv in factors]
	OUT.close()
	return()
	
def RunPEER(infile, covfile):
	print '1'
	import peer
	print '2'
	import scipy as SP
	print '3'
	expr = SP.loadtxt('%s/%s' % (cwd, '.'.join(infile.split('.')[:-1]) + '.NoHead.trans.tab'), delimiter='\t') # expr = SP.loadtxt('/home/vmason/LinuxShare/Programs/PEER/dataIN/DESeqNorm.untreated.NCBI.wVSD.2.tx.sort.trans.tab', delimiter='\t')
	expr.shape
	covs = SP.loadtxt('%s/%s' % (cwd, '.'.join(covfile.split('.')[:-1]) + '.NoHead.tab'), delimiter='\t') # covs = SP.loadtxt('/home/vmason/LinuxShare/Programs/PEER/dataIN/RAO_Covariates.no55.94.NoCTLRAO.NoMiss.trans.tab', delimiter='\t')
	covs.shape
	model = peer.PEER()
	model.setPhenoMean(expr)
	model.getPhenoMean().shape
	model.setCovariates(covs)
	model.getCovariates()
	model.setNk(10)
	model.getNk()
	#model.setAdd_mean(True) # the tutorial says that adding this value is usually a good idea
	#model.setTolerance(0.001) # default = 0.001
	model.getTolerance()
	#model.setVarTolerance(0.00001) # default = 0.00001
	model.getVarTolerance()
	model.setNmax_iterations(1000) # default = 1000

	model.update()
	
	
	#include covariates
	#covs = SP.loadtxt('/home/vmason/LinuxShare/Programs/PEER/dataIN/covariates.tab', delimiter='\t')
	#model.setCovariates(covs)
	
	#Model Parameters
	#model.setNmax_iterations(100)
	#model.setTolerance(0.01)
	#model.setVarTolerance(0.0001)
	#model.setPriorAlpha(0.001,0.1)
	#model.setPriorEps(0.1,10.)
	
	#In general you can keep the bound tolerance fairly high, but should keep the variation tolerance quite low compared to the variance of the expression matrix. If unsure, use the default values (bound=0.001, variance=0.00001).
	#PEER uses uninformative priors on weight precision and noise precision by default (Alpha a = 0.001, Alpha b = 0.1, Eps a = 0.1, Eps b = 10)
	
	
	factors = model.getX()
	weights = model.getW()
	precision = model.getAlpha()
	residuals = model.getResiduals()
	return(residuals, precision, factors, weights)

def FormatInputNoTranspose(f):
	o = []
	FILE = open("%s/%s" % (cwd, f), "r")
	header = FILE.readline()
	lines = FILE.readlines()
	FILE.close()
	
	names = [ l.strip().split(",")[0] for l in lines ]
	o = [ l.strip().split(",")[1:] for l in lines ]
	
	OUT = open("%s/%s" %(cwd, '.'.join(f.split('.')[:-1]) + '.NoHead.tab'), "w")
	OUT.write( '\n'.join( [ "\t".join(l) for l in o ] ) )
	OUT.close()
	
	return(header, names)
	
def FormatInput(f):
	o = []
	FILE = open("%s/%s" % (cwd, f), "r")
	header = FILE.readline()
	lines = FILE.readlines()
	FILE.close()
	
	names = [ l.strip().split("\t")[0] for l in lines ]
	o = [ l.strip().split("\t")[1:] for l in lines ]
	ot = map(list, zip(*o)) # requires all rows to be the same length or else the longer ones are removed
	
	OUT = open("%s/%s" %(cwd, '.'.join(f.split('.')[:-1]) + '.NoHead.trans.tab'), "w")
	OUT.write( '\n'.join( [ "\t".join(l) for l in ot ] ) )
	OUT.close()
	
	return(header, names)

def main():
	print 'Removing headers of gene expression matrix and transposing for Input into PEER'
	header, genenames = FormatInput(infile)
	print 'Removing headers of covariate file for Input into PEER'
	headercov, rownames = FormatInputNoTranspose(covfile)
	print 'Running PEER'
	residuals, precision, factors, weights = RunPEER(infile, covfile)
	print 'Writing out Residuals'
	WriteOUT(residuals, header, genenames)
	print 'Writing out Factors'
	WriteOUTFactors(factors, header)
	print 'Plot Factor Variance'
	plot0 = Plot_Alpha(precision)
	print 'plot Precision'
	plot1 = PlotGraph(precision)
	print 'plot Residuals'
	plot2 = PlotGraph(residuals)
	print 'plot Factors'
	plot3 = PlotGraph(factors)
	print 'plot Weights'
	plot4 = PlotGraph(weights)
	print 'save plots to: %s' % (pdffile)
	pp = PdfPages('%s/%s' % (cwd, pdffile))
	pp.savefig(plot0) #factor variance
	pp.savefig(plot1)
	pp.savefig(plot2)
	pp.savefig(plot3)
	pp.savefig(plot4)
	pp.close()
main()
