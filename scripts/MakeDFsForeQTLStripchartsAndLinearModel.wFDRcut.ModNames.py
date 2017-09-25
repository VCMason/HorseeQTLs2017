eqtls = 'D:\\LinuxShare\\Programs\\MatrixEQTL\\dataOUT\\MCK1_MAF0.05_VSD_no559443\\Cis.Linear.DESeqVSD.MCK1.c1e2t1e5.TagSNPs.NA.NoErrCov.PEERResid.NoNA.MAF0.05.best.tx.out'
covariates = 'D:\\LinuxShare\\Programs\\MatrixEQTL\\FormatDataIN\\Mock_wMAF0.05_VSD_no559443\\RAO_database.VCM.v4.MCK1.no559443.MatrixeQTL.csv'
geneexpression = 'D:\\LinuxShare\\Programs\\MatrixEQTL\\FormatDataIN\\Mock_wMAF0.05_VSD_no559443\\residuals.MCK1.VSD.MatrixeQTL.sort.tx.tab'
snps = 'D:\\LinuxShare\\Programs\\MatrixEQTL\\FormatDataIN\\Mock_wMAF0.05_VSD_no559443\\82_RAO_HD.2M.MCK1.trim2eQTLIndivs.ne1e3.MateQTL.TagSNP.tab'
#snpsloc = 'F:\LinuxShare\Programs\MatrixEQTL\FormatDataIN\HDE9_wMAF0.05_VSD\HDE9_RAO_HD.SNP.v6.sort.MAF0.05.NoNA.tab'
#geneloc = 'F:\LinuxShare\Programs\MatrixEQTL\FormatDataIN\HDE9_wMAF0.05_VSD\HDE9_RAO_HD.SNPLoc.v6.MatEQTL.trim.MAF0.05.tab'
fdrcut = 0.05

def MakePlot():
	blah = 'tahois'
	
def main():
	deqtls = {}
	names = []
	EQTLS = open(eqtls, 'r')
	line = EQTLS.readline() # skip header line
	line = EQTLS.readline()
	while line:
		s = line.strip().split('\t')
		if float(s[5]) < fdrcut:
			names.append('|'.join(s[:2]))
			#print '|'.join(s[:2])
			deqtls['|'.join(s[:2])] = s # key is snp name | gene name, value is list of line
		line = EQTLS.readline()
	EQTLS.close()
	
	covars = {}
	covarnames = []
	with open(covariates, 'r') as COVAR:
		lines = [line.strip().split(',') for line in COVAR]
	for x in lines: # x in list of list of lines from the covariate file 
		########## This is specific to my data ##########
		#modify individual names here
		if x[0] == 'BamFile':
			temp = []
			for i in x[1:]: # i == individual name
				if i[:2] == 'Un':
					n = '.'.join(i.split('_')[:2])
				else:
					n = i.split('_')[0]
				temp.append(n)
				#print n
			y = [x[0]] + temp
		else:
			y = x
		covars[x[0]] = list(y) # key is covariate name (individual) name
		covarnames.append(x[0])
		########## End Non-global code ##########
		########## Global code is below ##########
		#covars[x[0]] = list(x) # key is covariate name (individual) name
		#covarnames.append(x[0])
	
	genexp = {}
	with open(geneexpression, 'r') as GENEXP:
		lines = [line.strip().split('\t') for line in GENEXP]
	for x in lines:
		genexp[x[0]] = list(x) # key is gene name
		
	dsnp = {}
	with open(snps, 'r') as SNP:
		lines = [line.strip().split('\t') for line in SNP]
	for x in lines:
		dsnp[x[0]] = list(x) # key is snp name
	
	for snpgene in names: # deqtls.keys() # for each significant eqtl (snp|gene pair) export a file with columns: indiv, genotype, gene expression, and all covariates
		snp, gene = snpgene.split('|')
		output = []
		for i in range(len(covars[covarnames[0]])): # covarnames[0] == 'BamFile'
			covs = [ covars[c][i] for c in covarnames[1:] ]	# make list of all covariates (excluding first line)	
			output.append([covars[covarnames[0]][i], dsnp[snp][i], genexp[gene][i]] + covs) # add [ individual name, snp, gene exprssion, and all covariates ] to output in list for each individual (inside the big list output)
			#output.append([covars[covarnames[0]][i], dsnp[snp][i], genexp[gene][i]] + deqtls[snpgene][2:] + covs) # add [ individual name, snp, gene exprssion, beta, t-stat, p-value, FDR, and all covariates ] to output in list for each individual (inside the big list output)
		OUT = open(gene+'_'+snp+'.tab', 'w')
		OUT.write('\n'.join([ '\t'.join(x) for x in output]))
		OUT.close()
	
	
main()