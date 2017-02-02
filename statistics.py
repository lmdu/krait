#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PySide.QtSql import *

from db import *

class Statistics:
	meta_table = MetaTable()

class SequenceStatistics(Statistics):
	def __init__(self, unit='Mb', letter='ATCG'):
		self.table = FastaTable()
		self.unit = unit
		self.letter = letter
		self._bases = {'A':0,'G':0,'C':0,'T':0}
		self._total_sequences = 0
		self._total_bases = 0

	def _calc_bases(self):
		for fasta_file in self.table.fetchAll():
			fastas = pyfaidx.Fasta(fasta_file, sequence_always_upper=True)
			for fasta in fastas:
				seq = fasta[:].seq
				self._total_bases += len(seq)
				for b in self._bases:
					self._bases[b] += seq.count(b)

		for b in self._bases:
			self.meta_table.insert({'name': b, 'value': self._bases[b]})

		self.meta_table.insert({'name': 'total_bases', 'value': self._total_bases})

		return self._bases

	def _db_bases(self):
		for data in self.meta_table.query("SELECT * FROM meta WHERE name IN ('A', 'T', 'G', 'C')"):
			self._bases[data.name] = data.value

	@property
	def bases(self):
		if not self._bases:
			self._db_bases() or self._calc_bases()

		return self._bases

	def getSequenceCount(self):
		if not self._total_sequences:
			seq_table = SequenceTable()
			self._total_sequences = seq_table.get("SELECT COUNT(1) FROM sequence LIMIT 1")

		return self._total_sequences


	def getSize(self):
		'''
		get sequence total length A+G+T+C+N
		'''
		if not self._total_bases:
			self._total_bases = self.meta_table.getMeta('total_bases')

		if not self._total_bases:
			self._calc_bases()

		return self._total_bases

	def getNs(self):
		return self._total_bases - self.getValidSize()

	def getValidSize(self):
		'''
		get sequence valid total length A+G+T+C
		'''
		return sum(self.bases.values())

	def getGCContent(self):
		gc = self.bases['G'] + self.bases['C']
		return round(gc/float(self.getValidSize()), 2)

	def getScaleSize(self):
		scales = {'Mb': 1000000.0, 'Kb': 1000.0}
		if self.letter == 'ATGC':
			total = self.getValidSize()

		elif self.letter == 'ATGCN':
			total = self.getSize()

		return total/scales[self.unit]

	def getRelativeAbundance(self, counts):
		return round(counts/self.getScaleSize(), 2)
		
	def getRelativeDensity(self, lengths):
		return round(lengths/self.getScaleSize(), 2)

class MicrosatelliteStatistics(Statistics):
	def __init__(self):
		self.ssr_table = MicrosatelliteTable()

	@property
	def counts(self):
		return self.ssr_table.recordCounts()

	@property
	def lengths(self):
		return self.ssr_table.get("SELECT SUM(length) FROM ssr")

	def motifLength(self):
		sql = "SELECT length(motif) AS type, COUNT(1) AS count FROM ssr GROUP BY type"
		return {row.type: row.count for row in self.ssr_table.query(sql)}

	def motifType(self):
		sql = "SELECT standard, COUNT(1) AS count FROM ssr GROUP BY standard"
		return {row.standard: row.count for row in self.ssr_table.query(sql)}

	def motifRepeat(self):
		sql = "SELECT repeat, COUNT(1) AS count FROM ssr GROUP BY repeat"
		return {row.repeat: row.count for row in self.ssr_table.query(sql)}

	def report(self):
		pass

		



