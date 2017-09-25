import os
import glob
import subprocess
from datetime import datetime
startTime = datetime.now()

min_maf = '0.05'
window_len = '100'
min_r2 = '0.99'
max_len = '3'
merge_window_len = '100'
max_covered_times = '0'
mem_size = '0'
max_tagSNP_num = '0'

#outextsample = '_sample.txt' 
outextmaf = '.maf'
#outextmatrix = '.matrix'

def main():
	cwd = os.getcwd()
	files = glob.glob(cwd + "\\*" + outextmaf)
	#gens = glob.glob(cwd + "\\chr.+" + outextmatrix)
	for f in files:

		cmd = ' '.join([ 'FastTagger', '.'.join(f.split('\\')[-1].split('.')[:-1]), min_maf, window_len, min_r2, max_len, merge_window_len, max_covered_times, mem_size, max_tagSNP_num, '.'.join(f.split('\\')[-1].split('.')[:-1]) + '.out' ]) # Example: FastTagger ENm010 0.05  100 0.95  3  100  0  0  0  temp
		print cmd
		os.system(cmd)
		
	print datetime.now() - startTime
		
main()
