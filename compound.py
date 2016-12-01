#!/usr/bin/env python
import os
import sys
from db import SSRDB

def distance(ssr1, ssr2):
	return ssr2.start - ssr1.end - 1
	
def joinCompoundSSRs(cSSRs):
	chrom = cSSRs[0].chrom
	start = str(cSSRs[0].start)
	end = str(cSSRs[-1].end)
	complexity = str(len(cSSRs))

	motif = []
	cmotif = []
	seq = []
	cssrs = []

	for idx, ssr in enumerate(cSSRs):
		cmotif.append(ssr.cmotif)
		motif.append(ssr.motif)
		cssrs.append(ssr.id)
		if 0 < idx < len(cSSRs):
			d = distance(cSSRs[idx-1], ssr)
			if d > 0:
				seq.append("(N)%s" % d)
		seq.append("(%s)%s" % (ssr.motif, ssr.repeat))

	motif = "-".join(motif)
	cmotif = "-".join(cmotif)
	cssrs = ",".join(map(str,cssrs))
	seq = "-".join(seq)
	return (chrom, motif, cmotif, start, end, complexity, cssrs, seq)

def CompoundSSREngine(db, dMax):
	sql = "SELECT * FROM ssr ORDER BY chrom,start"
	ssrs = db.iter(sql)
	cSSRs = [ssrs.next()]

	for ssr in ssrs:
		d = distance(cSSRs[-1], ssr)
		if ssr.chrom == cSSRs[-1].chrom and d <= dMax:
			cSSRs.append(ssr)
		else:
			if len(cSSRs) > 1:
				yield joinCompoundSSRs(cSSRs)
			cSSRs = [ssr]
	else:
		if len(cSSRs) > 1:
			yield joinCompoundSSRs(cSSRs)


if __name__ == '__main__':
	db = SSRDB(sys.argv[1])
	for cSSR in CompoundSSREngine(db, 10):
		print "\t".join(cSSR)