#!/usr/bin/env python
from __future__ import division
import os
import re
import sys
import time
import pyfaidx
from motif import normalize
from utils import Data

__all__ = [
	"MicrosatelliteDetector",
	"CompoundDetector",
	"SatelliteDetector"
]

class MicrosatelliteDetector:
	'''
	Detect microsatellites from input fasta files with rules
	@para fasta str, input fasta file
	@para rules dict, minimal repeats for each ssr types
	'''
	def __init__(self, fasta, rules={1:12, 2:7, 3:5, 4:4, 5:4, 6:4}):
		self.fasta_file = fasta
		self.rules = rules

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
		return self.search()

	def _createFastaIndex(self):
		'''
		create fasta index using pyfaidx module
		'''
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
		#sequence processed progressing
		p1 = self.processed_sequences/self.total_sequences

		#sequence location progressing
		p2 = self.current_sequence_location/self.current_sequence_length
		return p1 * p2

	def search(self):
		'''
		search perfect microsatellites from the given fasta sequence file
		'''
		for fasta in self.fastas:
			self.processed_sequences += 1
			sequence = fasta[:].seq
			self.current_sequence_length = len(sequence)

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
				standard = normalize(motif)

				#get current sequence location
				self.current_sequence_location = stop

				yield Data(
					mid = None,
					sequence = fasta.name,
					start = start,
					stop = stop,
					motif = motif,
					standard = standard,
					type = mlen,
					repeat = repeat,
					length = length
				)


class CompoundDetector:
	'''
	detect compound microsatellites from perfect microsatellite results
	@para ssrs iterator, all SSRs
	@para dmax int, maximal allowed distance between two SSRs
	'''
	def __init__(self, ssrs, dmax=10):
		self.ssrs = ssrs
		self.dmax = dmax

	def __iter__(self):
		return self.search()

	def distance(self, ssr1, ssr2):
		return ssr2.start - ssr1.stop - 1

	def compound(self):
		motifs = []
		#standardized motifs
		standards = []
		#IDs of perfect SSRs
		components = []
		#structure of the compound SSRs
		structure = []
		#length of the total cSSRs
		length = 0

		for idx, ssr in enumerate(self.cssrs):
			standards.append(ssr.standard)
			motifs.append(ssr.motif)
			components.append(ssr.mid)
			length += ssr.length
			if 0 < idx < len(self.cssrs):
				d = self.distance(self.cssrs[idx-1], ssr)
				if d > 0:
					structure.append("(N)%s" % d)
			structure.append("(%s)%s" % (ssr.motif, ssr.repeat))
		
		return Data(
			cid = None,
			sequence = self.sequence,
			start = self.cssrs[0].start,
			stop = self.cssrs[-1].stop,
			motif = "-".join(motifs),
			standard = "-".join(standards),
			complexity = len(self.cssrs),
			length = length,
			component = ",".join(map(str, components)),
			structure = "-".join(structure)
		)

	def new(self, ssr):
		self.cssrs = [ssr]
		self.sequence = ssr.sequence

	def search(self):
		self.new(self.ssrs.next())
		for ssr in self.ssrs:
			d = self.distance(self.cssrs[-1], ssr)
			if ssr.sequence == self.sequence and d <= self.dmax:
				self.cssrs.append(ssr)
			else:
				if len(self.cssrs) > 1:
					yield self.compound()
				self.new(ssr)
		
		if len(cSSRs) > 1:
			yield self.compound()

class SatelliteDetector:
	'''
	Detect satellite DNA with any length of motif with any repeats
	@para motif tuple, minimal motif length and maximal motif length
	@para repeat int, the minimal tande repeats allowed
	'''
	def __init__(self, fasta, motifs=(7, 30), repeats=2):
		self.fasta_file = fasta
		self.motifs = motifs
		self.repeats = repeats

		#information of current processing fasta sequences
		self.total_sequences = 0
		self.processed_sequences = 0
		self.current_sequence_length = 0
		self.current_sequence_location = 0

		self._createSearchPattern()
		self._createFastaIndex()

	def __iter__(self):
		return self.search()

	def _createSearchPattern(self):
		'''
		create microsatellite search regular expression pattern
		'''
		self.pattern = re.compile(r'([ATGC]{1,%s}?)\1{%s,}' % (self.motifs[1], self.repeats-1))

	def _createFastaIndex(self):
		if not os.path.exists(self.fasta_file):
			raise Exception('Fasta file %s is not exists' % self.fasta_file)
		
		self.fastas = pyfaidx.Fasta(self.fasta_file, sequence_always_upper=True)
		
		#total sequences in fasta file
		self.total_sequences = len(self.fastas.keys())

	@property
	def progress(self):
		'''
		Get current process and return float 0-1
		'''
		p1 = self.processed_sequences/self.total_sequences
		p2 = self.current_sequence_location/self.current_sequence_length
		return p1 * p2

	def search(self):
		'''
		search perfect microsatellites from the given fasta file
		'''
		for fasta in self.fastas:
			self.processed_sequences += 1
			sequence = fasta[:].seq
			self.current_sequence_length = len(seq)
			
			for tandem in self.pattern.finditer(sequence):
				ssr, motif = tandem.group(0,1)
				#get motif length
				mlen = len(motif)

				if mlen < self.motifs[0] or mlen > self.motifs[1]:
					continue

				#get length of ssr sequence
				length = len(ssr)
				
				#get repeats of the ssr
				repeat = length/mlen

				#get the start and end location of ssr
				start, stop = tandem.span()
				start += 1

				#get current sequence location
				self.current_sequence_location = stop

				yield Data(
					sid = None,
					sequence = name,
					start = start,
					stop = stop,
					motif = motif,
					type = mlen,
					repeat = repeat,
					length = length
				)


if __name__ == '__main__':
	start_time = time.time()
	for ssr in MicrosatelliteDetector(sys.argv[1]):
		pass
	end_time = time.time()
	print round(end_time - start_time, 2)
	
