#!/usr/bin/env python
import fasta
import tandem
import time
seqs = fasta.GzipFasta('dr.fa')

def test():
	with seqs:
		for name, seq in seqs:
			ssr = tandem.search_ssr(seq, (12,7,5,4,4,4))
			print name
			#ssr = tandem.search_issr(seq, 3,8,2,1,2,12,2000)

test()
