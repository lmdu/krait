#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import json

from libs import *
from db import *
from utils import *
from config import *

class Statistics(object):
	_db = None
	_total_bases = 0
	_total_sequences = 0
	_bases = {}
	unit = None
	letter = None
	def __init__(self, unit='Mb', letter='ATGC'):
		self.unit = unit
		self.letter = letter

	@property
	def db(self):
		if self._db is None:
			self._db = Database()
		return self._db

	def count_bases(self):
		for fasta_file in self.db.get_column("SELECT path FROM fasta"):
			fastas = fasta.GzipFasta(fasta_file)
			for name in fastas.keys:
				sequence = fastas[name]
				self._total_bases += len(sequence)
				self._total_sequences += 1
				for b in ['A','T','G','C']:
					self._bases[b] = self._bases.get(b, 0) + sequence.count(b)

		for b in self.bases:
			self.db.set_option('statis_base_%s' % b, self._bases[b])

		self.db.set_option('statis_base_total', self._total_bases)
		self.db.set_option('statis_seqs_count', self._total_sequences)

	def get_bases(self):
		for b in ['A', 'T', 'G', 'C']:
			self._bases[b] = self.db.get_option('statis_base_%s' % b)

	@property
	def bases(self):
		if not self._bases:
			self.get_bases()
		
		return self._bases

	@property
	def seqcount(self):
		if not self._total_sequences:
			seq_counts = self.db.get_option('statis_seqs_count')
			if seq_counts:
				self._total_sequences = seq_counts
			else:
				self.count_bases()

		return self._total_sequences

	@property
	def size(self):
		if not self._total_bases:
			self._total_bases = self.db.get_option('statis_base_total')

		return self._total_bases

	@property
	def validsize(self):
		return sum(self.bases.values())

	@property
	def ns(self):
		return self.size - self.validsize

	@property
	def transize(self):
		'''
		get transformed size bp to Mb or bp to Kb
		'''
		scales = {'Mb': 1000000, 'Kb': 1000}
		if self.letter == 'ATGC':
			total = self.validsize

		elif self.letter == 'ATGCN':
			total = self.size

		return total/scales[self.unit]

	@property
	def gc(self):
		'''
		get the GC content of the fastas
		'''
		gc = self.bases['G'] + self.bases['C']
		return round(gc/self.validsize*100, 2)

	def ra(self, counts):
		'''
		calculate the relative abundance of given counts
		'''
		return round(counts/self.transize, 2)
		
	def rd(self, lengths):
		'''
		calculate the relative density of given lengths
		'''
		return round(lengths/self.transize, 2)

	def type(self, index):
		types = {1: 'Mono', 2: 'Di', 3: 'Tri', 4: 'Tetra', 5: 'Penta', 6: 'Hexa'}
		return types[index]

	def region(self, cat):
		total_counts = self.db.get_one("SELECT COUNT(1) FROM %s LIMIT 1" % cat)
		sql = "SELECT feature,COUNT(1) AS count FROM location WHERE category='%s' GROUP BY feature"
		rows = []
		feat_counts = 0
		for row in self.db.query(sql % cat):
			rows.append((row.feature, row.count))
			feat_counts += row.count
		
		if rows:
			rows.append(('Intergenic', total_counts-feat_counts))
		
		return rows

	def results(self):
		r = Data()
		r.seqcount = self.seqcount
		r.size = self.size
		r.validsize = self.validsize
		r.ns = self.ns, round(self.ns/self.size*100, 2)
		r.gc = self.gc
		r.unit = self.unit
		return r


class SSRStatistics(Statistics):
	_ssr_counts = 0
	_ssr_length = 0
	
	def __init__(self):
		super(SSRStatistics, self).__init__()

	@property
	def count(self):
		'''
		get total perfect SSR counts and frequency
		'''
		if not self._ssr_counts:
			self._ssr_counts = self.db.get_one("SELECT COUNT(1) FROM ssr LIMIT 1")
		
		return self._ssr_counts

	@property
	def length(self):
		if not self._ssr_length:
			self._ssr_length = self.db.get_one("SELECT SUM(length) FROM ssr LIMIT 1")
		
		return self._ssr_length

	@property
	def frequency(self):
		return self.ra(self.count)

	@property
	def density(self):
		return self.rd(self.length)
		
	def motifTypeStatis(self):
		sql = "SELECT type, SUM(length) AS length, COUNT(1) AS count FROM ssr GROUP BY type ORDER BY type"
		rows = [('Type', 'Counts', 'Length (bp)', 'Percent (%)', 'Average Length (bp)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.db.query(sql):
			percent = round(row.count/self.count*100, 2)
			average = round(row.length/row.count, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((self.type(row.type), row.count, row.length, percent, average, frequency, density))
		return rows

	def motifCategoryStatis(self):
		sql = "SELECT standard, SUM(length) AS length, COUNT(1) AS count FROM ssr GROUP BY standard ORDER BY type,standard"
		rows = [('Motif', 'Counts', 'Length (bp)', 'Percent (%)', 'Average Length (bp)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.db.query(sql):
			percent = round(row.count/self.count*100, 2)
			average = round(row.length/row.count, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((row.standard, row.count, row.length, percent, average, frequency, density))
		return rows

	def motifRepeatStatis(self):
		rows = []
		for i in range(1,7):
			sql = "SELECT repeat FROM ssr WHERE type=%s" % i
			r = self.db.get_column(sql)
			if r: rows.append(r)
		return rows
		
	def SSRLengthStatis(self):
		rows = []
		for i in range(1,7):
			sql = "SELECT length FROM ssr WHERE type=%s" % i
			r = self.db.get_column(sql)
			if r: rows.append(r)
		return rows

	def results(self):
		r = Data()
		r.count = self.count
		r.length = self.length
		r.frequency = self.frequency
		r.density = self.density
		r.type = self.motifTypeStatis()
		r.category = self.motifCategoryStatis()
		r.repeat = self.motifRepeatStatis()
		r.ssrlen = self.SSRLengthStatis()
		r.region = self.region('ssr')
		return r

class CSSRStatistics(Statistics):
	_cssr_counts = 0
	_cm_counts = 0
	_cssr_length = 0
	def __init__(self):
		super(CSSRStatistics, self).__init__()

	@property
	def cssr_count(self):
		if not self._cssr_counts:
			self._cssr_counts = self.db.get_one("SELECT SUM(complexity) FROM cssr LIMIT 1")

		return self._cssr_counts

	@property
	def cm_count(self):
		if not self._cm_counts:
			self._cm_counts = self.db.get_one("SELECT COUNT(1) FROM cssr LIMIT 1")

		return self._cm_counts

	@property
	def length(self):
		if not self._cssr_length:
			self._cssr_length = self.db.get_one("SELECT SUM(length) FROM cssr LIMIT 1")

		return self._cssr_length

	@property
	def frequency(self):
		return self.ra(self.cm_count)

	@property
	def density(self):
		return self.rd(self.length)

	def CSSRComplexityStatis(self):
		sql = "SELECT complexity,COUNT(1) FROM cssr GROUP BY complexity ORDER BY complexity"
		rows = [('Complexity', 'Counts')]
		for row in self.db.query(sql):
			rows.append(row.values)
		return rows

	def CSSRLengthStatis(self):
		sql = "SELECT length,COUNT(1) FROM cssr GROUP BY length ORDER BY length"
		rows = [('Length', 'Counts')]
		for row in self.db.query(sql):
			rows.append(row.values)
		return rows

	def CSSRGapStatis(self):
		sql = "SELECT gap,COUNT(1) FROM cssr GROUP BY gap ORDER BY gap"
		rows = [('Gap length', 'Counts')]
		for row in self.db.query(sql):
			rows.append(row.values)
		return rows

	def results(self):
		r = Data()
		r.cssr_count = self.cssr_count
		r.cm_count = self.cm_count
		r.length = self.length
		r.frequency = self.frequency
		r.density = self.density
		r.complexity = self.CSSRComplexityStatis()
		r.cssrlen = self.CSSRLengthStatis()
		r.gap = self.CSSRGapStatis()
		r.region = self.region('cssr')
		return r

class ISSRStatistics(Statistics):
	_issr_counts = 0
	_issr_length = 0
	def __init__(self):
		super(ISSRStatistics, self).__init__()

	@property
	def count(self):
		if not self._issr_counts:
			self._issr_counts = self.db.get_one("SELECT COUNT(1) FROM issr LIMIT 1")

		return self._issr_counts

	@property
	def length(self):
		if not self._issr_length:
			self._issr_length = self.db.get_one("SELECT SUM(length) FROM issr LIMIT 1")

		return self._issr_length

	@property
	def frequency(self):
		return self.ra(self.count)

	@property
	def density(self):
		return self.rd(self.length)

	def motifTypeStatis(self):
		sql = "SELECT type, SUM(length) AS length, COUNT(1) AS count FROM issr GROUP BY type ORDER BY type"
		rows = [('Type', 'Counts', 'Length (bp)', 'Percent (%)', 'Average Length (bp)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.db.query(sql):
			percent = round(row.count/self.count*100, 2)
			average = round(row.length/row.count, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((self.type(row.type), row.count, row.length, percent, average, frequency, density))
		return rows

	def motifCategoryStatis(self):
		sql = "SELECT standard, SUM(length) AS length, COUNT(1) AS count FROM issr GROUP BY standard ORDER BY type,standard"
		rows = [('Motif', 'Counts', 'Length (bp)', 'Percent (%)', 'Average Length (bp)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.db.query(sql):
			percent = round(row.count/self.count*100, 2)
			average = round(row.length/row.count, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((row.standard, row.count, row.length, percent, average, frequency, density))
		return rows

	def ISSRScoreStatis(self):
		rows = []
		for i in range(1,7):
			sql = "SELECT score FROM issr WHERE type=%s" % i
			r = self.db.get_column(sql)
			if r: rows.append(r)
		return rows

	def ISSRLengthStatis(self):
		rows = []
		for i in range(1,7):
			sql = "SELECT length FROM issr WHERE type=%s" % i
			r = self.db.get_column(sql)
			if r: rows.append(r)
		return rows

	def results(self):
		r = Data()
		r.count = self.count
		r.length = self.length
		r.frequency = self.frequency
		r.density = self.density
		r.type = self.motifTypeStatis()
		r.category = self.motifCategoryStatis()
		r.score = self.ISSRScoreStatis()
		r.issrlen = self.ISSRLengthStatis()
		r.region = self.region('issr')
		return r

class VNTRStatistics(Statistics):
	_vntr_counts = 0
	_vntr_length = 0
	def __init__(self):
		super(VNTRStatistics, self).__init__()

	@property
	def count(self):
		if not self._vntr_counts:
			self._vntr_counts = self.db.get_one("SELECT COUNT(1) FROM vntr LIMIT 1")
		return self._vntr_counts

	@property
	def length(self):
		if not self._vntr_length:
			self._vntr_length = self.db.get_one("SELECT SUM(length) FROM vntr LIMIT 1")
		return self._vntr_length

	@property
	def frequency(self):
		return self.ra(self.count)

	@property
	def density(self):
		return self.rd(self.length)

	def motifTypeStatis(self):
		sql = "SELECT type, COUNT(1) AS count FROM vntr GROUP BY type"
		rows = [row.values for row in self.db.query(sql)]
		return rows

	def motifRepeatStatis(self):
		sql = "SELECT repeat, COUNT(1) AS count FROM vntr GROUP BY repeat"
		rows = [row.values for row in self.db.query(sql)]
		return rows

	def VNTRLengthStatis(self):
		sql = "SELECT length, COUNT(1) AS count FROM vntr GROUP BY length"
		rows = [row.values for row in self.db.query(sql)]
		return rows

	def results(self):
		r = Data()
		r.count = self.count
		r.length = self.length
		r.frequency = self.frequency
		r.density = self.density
		r.type = self.motifTypeStatis()
		r.repeat = self.motifRepeatStatis()
		r.vntrlen = self.VNTRLengthStatis()
		r.region = self.region('vntr')
		return r


