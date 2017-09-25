import os
import glob
import pandas
from datetime import datetime
startTime = datetime.now()

span = 4

extsample = '_samples.txt' 
extmatrix = '.matrix'
extMAFSNPs = '.out.SNPlineno.txt' # contains number on each line, contains the number MINUS 1 of the SNPs in file ".maf" and ".matrix" that satisfies the min_maf threshold. (i.e. start counting from 0)
extsnpid = '.out.SNPid.txt' # contains SNPs names/info that passed MAF filter
exttag = '.out.tagSNP.txt' # contains number on each line, number refers to SNP on number-1 line of .snpid.txt (line number of SNPs that pass MAF threshold), these SNPs are TAG SNPs (i.e. start counting from 0)

def WriteOUT(f, output):
	OUT = open(f, 'wb')
	OUT.write(output)
	OUT.close()
	
def FormatForeQTLBMA(maf, gens012, samples, prefix):
	header = '\t'.join(samples) + '\n'
	out = '\n'.join([ '\t'.join([ i.strip().split('\t')[0] ] + j) for i, j in map(None, maf, gens012) ]) # i.strip().split('\t')[0] is the SNP name, j is a list of all genotypes for one SNP coded as 0,1,2
	WriteOUT(prefix + '.eQTLBMA.SNP.tab', header + out)
	outloc = '\n'.join([ '\t'.join([ prefix.split('.')[-1], i.strip().split('\t')[1], str(int(i.strip().split('\t')[1])+1), i.strip().split('\t')[0] ]) for i in maf ]) # chr, pos, pos+1, snpid
	WriteOUT(prefix + '.eQTLBMA.SNPLoc.tab', out)
	return(header, out, outloc)
	
def FormatForMatrixeQTL(maf, gens012, samples, prefix):
	header = '\t'.join(['snpid'] + samples)  + '\n'
	out = '\n'.join([ '\t'.join([ i.strip().split('\t')[0] ] + j) for i, j in map(None, maf, gens012) ]) # i.strip().split('\t')[0] is the SNP name, j is a list of all genotypes for one SNP coded as 0,1,2
	WriteOUT(prefix + '.MateQTL.SNP.tab', header + out)
	headerloc = 'snpid\tchr\tpos\n'
	outloc = '\n'.join([ '\t'.join([ i.strip().split('\t')[0], prefix.split('.')[-1], i.strip().split('\t')[1] ]) for i in maf ]) # snpid, chr, pos
	WriteOUT(prefix + '.MateQTL.SNPLoc.tab', headerloc + out)
	return(header, headerloc, out, outloc)

def Call012ForGenotypes(matrix): # collapses biallelic genotypes to 0 for homozygour reference allele, to 1 for heteroygous, and 2 for homozygous alternative allele
	gens012 = [ [ '0' if snp[i:i+span].strip() == '0 0' else '1' if snp[i:i+span].strip() == '0 1' else '1' if snp[i:i+span].strip() == '1 0' else '2' if snp[i:i+span].strip() == '1 1' else 'something is:%s: wrong' % (snp[i:i+span].strip()) for i in range(0,len(snp),span) ] for snp in matrix ] # list of list, outer list is list of all SNPs, inner list of 0,1,2 coded genotypes (one for each individual for that SNP)
	return(gens012)
	
def FilterMatrixByIndices(matrix, indices): # FilterMatrixForMAFSNPs # and/or FilterMatrixForTAGSNPs
	filteredmatrix = [ matrix[int(i)] for i in indices ]
	return(filteredmatrix)
	
def main():
	cwd = os.getcwd()
	matrices = glob.glob(cwd + "\\*" + extmatrix)
	
	allsnpMateQTL = ''
	allsnplocMateQTL = ''
	allsnpeQTLBMA = ''
	allsnploceQTLBMA = ''
	count = 0
	for fmatrix in matrices:
		prefix = fmatrix[:-len(extmatrix)]
		fsample, fMAFSNPs, fsnpid, ftag  = prefix + extsample, prefix + extMAFSNPs, prefix + extsnpid, prefix + exttag
		
		F = open(fmatrix, 'r')
		matrix = F.readlines()
		F.close()		
		F = open(fsample, 'r')
		samples = F.readlines()
		F.close()
		F = open(fMAFSNPs, 'r')
		MAFSNPIndices = F.readlines()
		F.close()
		F = open(fsnpid, 'r')
		MAFsnpids = F.readlines()
		F.close()
		F = open(ftag, 'r')
		tagSNPIndices = F.readlines()
		F.close()
		
		samples = [ s.strip() for s in samples ]
		MAFmatrix = FilterMatrixByIndices(matrix, MAFSNPIndices)
		TAGmatrix = FilterMatrixByIndices(MAFmatrix, tagSNPIndices)
		TAGSNPids = FilterMatrixByIndices(MAFsnpids, tagSNPIndices)
		gens012 = Call012ForGenotypes(TAGmatrix)
		MateQTLheader, MateQTLheaderloc, tempMateQTLsnp, tempMateQTLsnploc = FormatForMatrixeQTL(TAGSNPids, gens012, samples, prefix)
		eQTLBMAheader, tempeQTLBMAsnp, tempeQTLBMAsnploc = FormatForeQTLBMA(TAGSNPids, gens012, samples, prefix)
		if count > 0:
			allsnpMateQTL += tempMateQTLsnp + '\n'
			allsnplocMateQTL += tempMateQTLsnploc + '\n'
			allsnpeQTLBMA += tempeQTLBMAsnp + '\n'
			allsnploceQTLBMA += tempeQTLBMAsnploc + '\n'
		elif count == 0:
			allsnpMateQTL += MateQTLheader + tempMateQTLsnp + '\n'
			allsnplocMateQTL += MateQTLheaderloc + tempMateQTLsnploc + '\n'
			allsnpeQTLBMA += eQTLBMAheader + tempeQTLBMAsnp + '\n'
			allsnploceQTLBMA += tempeQTLBMAsnploc + '\n'
		
		count += 1

	WriteOUT('.'.join(prefix.split('.')[:-1]) + '.MateQTL.AllSNP.tab', allsnpMateQTL)
	WriteOUT('.'.join(prefix.split('.')[:-1]) + '.MateQTL.AllSNPLoc.tab', allsnplocMateQTL)
	WriteOUT('.'.join(prefix.split('.')[:-1]) + '.eQTLBMA.AllSNP.tab', allsnpeQTLBMA)
	WriteOUT('.'.join(prefix.split('.')[:-1]) + '.eQTLBMA.AllSNPLoc.tab', allsnploceQTLBMA)

	print datetime.now() - startTime
		
main()
