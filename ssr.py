#!/usr/bin/env python
import os
import sys
import time
import pyfaidx
from motif import MOTIFS
from optparse import OptionParser

import re

#Search all ssrs from sequence
#@para sequence string, fasta sequence.
#@para threshold dict, the minimal repeats for ssr.
#@return tuple, ssr information
def SSRSearchEngine(chrom, sequence, threshold):
	pattern = re.compile(r'([atgcATGC]{1,6}?)\1{%s,}' % (min(threshold.values()) - 1))
	ssrs = pattern.finditer(sequence)
	for ssr in ssrs:
		seq = ssr.group(0)
		length = len(seq)
		motif = ssr.group(1)
		motifLen = len(motif)
		repeat = length/motifLen
		if repeat < threshold[motifLen]: continue
		start = ssr.start(0) + 1
		end = ssr.end(0)
		smotif = MOTIFS[motif.upper()]
		yield (None, chrom, start, end, repeat, length, motif, smotif)

def main():
	threshold = {1:12, 2:7, 3:5, 4:4, 5:4, 6:4}
	start_time = time.time()
	genome = pyfaidx.Fasta('../chr1.fa')
	for chrom in genome:
		for ssr in SSRSearchEngine(chrom.name, str(chrom[:].seq), threshold):
			pass
	end_time = time.time()
	print round(end_time - start_time, 2)

if __name__ == '__main__':
	import cProfile
	cProfile.run("main()")
	
