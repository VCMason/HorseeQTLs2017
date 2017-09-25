'''
Much of this code was lifted from Damien Farrell: https://github.com/dmnfarrell
Author: Victor Mason
Input: Gene expression count matrix, individuals as columns rows are genes, and tab delimited.
Input: First header line should only have individual names, and no geneid header.
id1 id2 id3
MC1R 5 20 14
ASIP 10 13 8
'''

import numpy as np
import pandas as pd
import matplotlib

f = "D:/LinuxShare/Programs/Python/KS_Test/MCK1/DESeqNorm.NCBI.NoVSD.RSum1Trim.MCK1.tx.tab" # path to DESeq normalized=TRUE read counts.
outext = '.RowSum1trim.norm'
diff = 0.0001


def filterExprResults(df, rc):
	df = df.loc[df.mean(axis=1)>=rc,:] # keeps only rows (i.e. genes) with mean read count >= mean read count cutoff
	return(df)
	
def ks_test(df, ids): # no need for pairs of ids, i will create them below.
	''' 
		Computes KS statistic, returns d-statistic, and pvaue for each pair of individuals.
		df is a pandas dataframe of the gene expression file
		ids is a list of all possible pairs of individual names
	''' 
	def getMin(x):
		r=0
		x=list(x)
		m=x[0]
		#print x
		#print len(x)
		for i in range(len(x)):
			if x[i]<m and x[i+1]>x[i] and abs(m-x[i])>diff: # current D-statistic must be lower than starting D-statistic, next D-statistic must be large than current, current d-statistic must be less than 0.02 lower than starting D-statistic
				r=readcuts[i]
				break
		return r
	
	from scipy.stats import ks_2samp
	import pylab as plt
	
	readcuts = np.arange(0,300,2) # make array of cuttoff values
	result = []
	count = 0
	for i in ids:
		stats = []
		s1, s2 = i
		for rc in readcuts:
			d = filterExprResults(df, rc)
			x = d[s1]
			y = d[s2]
			val = ks_2samp(x, y) # val is D-statistic, and pvalue
			stats.append(val[0]) # grab only the D-statistic
			#print s1, s2, str(rc)
			#print val
		result.append(stats)
		count += 1
		if count % 100 == 0:
			print count
	result = pd.DataFrame(np.array(result).T)
	result=result.set_index(readcuts)
	outids = ["|".join(i) for i in ids]
	result.columns = outids
	result.to_csv('Dstats%s_%f.txt' % (outext, diff), index=True, header=True, sep='\t')
	tmins = result.apply(getMin) # sends one columns at a time, one column is one pair of samples with a D-statistic for each read count cutoff.
	#print 'These are the cutoff values:'
	#print (tmins)
	mean=tmins.mean()
	median=tmins.median()
	print 'Mean read cutoff value for all pairs of samples: %d' % (mean)
	print 'Median read cutoff value for all pairs of samples: %d' % (median)
	tmins.to_csv('AllMeanReadCutoffValues%s_%f.txt' % (outext, diff), index=True, header=True, sep='\t')
	
	std=tmins.std()
	
	result.plot(lw=0.5,colormap='Set1',legend=False)

	plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
	plt.xlabel('read count threshold')
	plt.ylabel('KS')
	plt.savefig('KS_test.AllSamples%s_%f.pdf' % (outext, diff))
	
	randcols = result.sample(10, axis=1)
	randcols.plot(lw=1,colormap='Set1',legend=False)

	# Make plot readable
	plt.axvspan(mean-std/2,mean+std/2,color='g',alpha=0.3)
	plt.xlabel('read count threshold')
	plt.ylabel('KS')
	plt.savefig('KS_test.5Samples%s_%f.pdf' % (outext, diff))
	
def CreateListOfAllIDPairs(names):
	'''
	Example:
	for i in itertools.combinations([1,2,3], 2):
		print i
	'''
	import itertools
	idpairs = []
	for i in itertools.combinations(names, 2):
		idpairs.append(i)
	return(idpairs)
	
def main():
	print 'Reading gene expression file'	
	df = pd.read_csv(f, sep='\t') # pandas dataframe of expression values
	
	names = list(df.columns.values)
	
	print 'Creating all pairs of sample IDs'
	ids = CreateListOfAllIDPairs(names)
	print 'Number of sample IDs; %d\n Number of pairs of samples: %d' % (len(names), len(ids))
	
	print 'Starting KS-Test iteration'
	ks_test(df, ids)
	
main()