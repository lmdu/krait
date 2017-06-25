#!/usr/bin/env python
import pyfaidx
import tandem

f = pyfaidx.Fasta('chr13.fa')
seq = f['chr13'][:].seq.upper()

issrs = tandem.search_issr(seq, 3, 8, 3, 50, 200)

for issr in issrs:
	s = seq[issr[2]-1:issr[3]]
	issr[-1] = round(issr[-1], 2)
	print "%s\t%s" % ("\t".join(map(str, issr)), s)
