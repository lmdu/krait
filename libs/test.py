#!/usr/bin/env python
import pyfaidx
import tandem

f = pyfaidx.Fasta('chr13.fa')
seq = f['chr13'][:].seq.upper()

issrs = tandem.search_issr(seq, 3, 8, 3, 12)

for issr in issrs:
	s = seq[issr[2]-1:issr[3]]
	print "%s\t%s" % ("\t".join(map(str, issr)), s)
