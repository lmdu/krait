#!/usr/bin/env python
import os
import re
import sys
import time
import pyfaidx
from motif import normalize

class SSR(dict):
	def __getattr__(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

	def __setattr__(self, name, val):
		self[name] = val


class SSRDetector:
	'''
	Detect microsatellites from input fasta files with rules
	@para fasta str, input fasta file
	@para rules dict, minimal repeats for each ssr types
	'''
	def __init__(self, fasta, rules=None):
		self.fasta_file = fasta
		self.rules = rules

		#default rules
		if rules is None:
			self.rules = {1:12, 2:7, 3:5, 4:4, 5:4, 6:4}

		#create ssr search pattern
		self._createSearchPattern()

	def __iter__(self):
		for name, seq in self.getSequences():
			for ssr in self.searchSSRsInSequence(name, seq):
				yield ssr

	def createFastaIndex(self):
		if not os.path.exists(self.fasta_file):
			raise Exception('Fasta file %s is not exists' % self.fasta_file)
		
		self.fasta = pyfaidx.Fasta(self.fasta_file)

	def _createSearchPattern(self):
		'''
		create microsatellite search regular expression pattern
		'''
		repeats = min(self.rules.values()) - 1
		self.pattern = re.compile(r'([ATGC]{1,6}?)\1{%s,}' % repeats)

	def getSequenceCount(self):
		return len(self.fasta.keys())
	
	def getSequences(self):
		for seq in self.fasta:
			yield (seq.name, seq[:].seq.upper())

	def searchSSRsInSequence(self, chrom, sequence):
		for tandem in self.pattern.finditer(sequence):
			ssr, motif = tandem.group(0,1)
			length = len(ssr)
			mlen = len(motif)
			repeat = length/mlen
			
			#ssr can not match the search rule
			if repeat < self.rules[mlen]:
				continue

			#get the start and end location of ssr
			start, end = tandem.span()
			start += 1

			#standardize the motif of ssr
			smotif = normalize(motif)

			yield SSR(
				ID = None,
				sequence = chrom,
				start = start,
				end = end, 
				repeat = repeat,
				length = length,
				motif = motif,
				smotif = smotif
			)

def join_compound_SSRs(cSSRs):
	sequence = cSSRs[0].sequence
	start = str(cSSRs[0].start)
	end = str(cSSRs[-1].end)
	complexity = str(len(cSSRs))

	motif = []
	smotif = []
	seq = []
	cssrs = []

	for idx, ssr in enumerate(cSSRs):
		smotif.append(ssr.smotif)
		motif.append(ssr.motif)
		cssrs.append(ssr.ID)
		if 0 < idx < len(cSSRs):
			d = distance(cSSRs[idx-1], ssr)
			if d > 0:
				seq.append("(N)%s" % d)
		seq.append("(%s)%s" % (ssr.motif, ssr.repeat))

	motif = "-".join(motif)
	smotif = "-".join(smotif)
	cssrs = ",".join(map(str,cssrs))
	seq = "-".join(seq)
	return SSR(
		sequence = sequence,
		start = start,
		end = end,
		motif = motif,
		smotif = smotif,
		complexity = complexity,
		cssrs = cssrs,
		compound = seq
	)

def distance(ssr1, ssr2):
	return ssr2.start - ssr1.end - 1	

class CSSRDetector:
	def __init__(self, ssrs, dmax=None):
		self.ssrs = ssrs
		self.dmax = dmax
		if dmax is None:
			self.dmax = 10

	def __iter__(self):
		return self.searchCompound()

	def searchCompound(self):
		cSSRs = [self.ssrs.next()]
		for ssr in self.ssrs:
			d = distance(cSSRs[-1], ssr)
			if ssr.sequence == cSSRs[-1].sequence and d <= self.dmax:
				cSSRs.append(ssr)
			else:
				if len(cSSRs) > 1:
					yield join_compound_SSRs(cSSRs)
				cSSRs = [ssr]

		else:
			if len(cSSRs) > 1:
				yield join_compound_SSRs(cSSRs)


	

if __name__ == '__main__':
	start_time = time.time()
	for ssr in MicrosatelliteDetector(sys.argv[1]):
		pass
	end_time = time.time()
	print round(end_time - start_time, 2)
	
