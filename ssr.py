#!/usr/bin/env python
from __future__ import division
import os
import re
import sys
import time
import pyfaidx
from motif import normalize
from utils import Data

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

		#information of current processing fasta sequences
		self.total_sequences = 0
		self.processed_sequences = 0
		self.current_sequence_length = 0
		self.current_sequence_location = 0

		#create ssr search pattern
		self._createSearchPattern()

		#create fasta faidx index
		self._createFastaIndex()

	def __iter__(self):
		for fasta in self.fastas:
			self.processed_sequences += 1
			seq = str(fasta)
			self.current_sequence_length = len(seq)
			for ssr in self.search(fasta.name, seq):
				yield ssr

	def _createFastaIndex(self):
		if not os.path.exists(self.fasta_file):
			raise Exception('Fasta file %s is not exists' % self.fasta_file)
		
		self.fastas = pyfaidx.Fasta(self.fasta_file, sequence_always_upper=True)
		
		#total sequences in fasta file
		self.total_sequences = len(self.fastas.keys())

	def _createSearchPattern(self):
		'''
		create microsatellite search regular expression pattern
		'''
		repeats = min(self.rules.values()) - 1
		self.pattern = re.compile(r'([ATGC]{1,6}?)\1{%s,}' % repeats)
	
	@property
	def progress(self):
		'''
		Get current process and return float 0-1
		'''
		p1 = self.processed_sequences/self.total_sequences
		p2 = self.current_sequence_location/self.current_sequence_length
		return p1 * p2

	def search(self, name, sequence):
		'''
		search perfect microsatellites from the given sequence
		@para name str, sequence name
		@para sequence str, DNA sequence with ATGC
		@return dict, one ssr with base information
		'''
		for tandem in self.pattern.finditer(sequence):
			ssr, motif = tandem.group(0,1)
			#get length of ssr sequence
			length = len(ssr)
			#get motif length
			mlen = len(motif)
			#get repeats of the ssr
			repeat = length/mlen
			
			#ssr can not match the search rule
			if repeat < self.rules[mlen]:
				continue

			#get the start and end location of ssr
			start, stop = tandem.span()
			start += 1

			#standardize the motif of ssr
			smotif = normalize(motif)

			#get current sequence location
			self.current_sequence_location = stop

			yield Data(
				ID = None,
				sequence = name,
				start = start,
				stop = stop,
				motif = motif,
				smotif = smotif,
				mlength = mlen,
				repeat = repeat,
				length = length
			)

def join_compound_SSRs(cSSRs):
	sequence = cSSRs[0].sequence
	start = cSSRs[0].start
	stop = cSSRs[-1].stop
	complexity = len(cSSRs)

	motif = []
	smotif = []
	seq = []
	cssrs = []
	size = 0

	for idx, ssr in enumerate(cSSRs):
		smotif.append(ssr.smotif)
		motif.append(ssr.motif)
		cssrs.append(ssr.ID)
		size += ssr.length
		if 0 < idx < len(cSSRs):
			d = distance(cSSRs[idx-1], ssr)
			if d > 0:
				seq.append("(N)%s" % d)
		seq.append("(%s)%s" % (ssr.motif, ssr.repeat))

	motif = "-".join(motif)
	smotif = "-".join(smotif)
	cssrs = ",".join(map(str,cssrs))
	seq = "-".join(seq)
	return Data(
		ID = None,
		sequence = sequence,
		start = start,
		stop = stop,
		motif = motif,
		smotif = smotif,
		complexity = complexity,
		length = size,
		cssrs = cssrs,
		compound = seq
	)

def distance(ssr1, ssr2):
	return ssr2.start - ssr1.stop - 1	

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
	
