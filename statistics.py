#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import json
import pygal
import pyfaidx

from pygal.style import DefaultStyle
from PySide.QtCore import QDir

from db import *
from utils import Data

class Plots:
	def __init__(self, name, data):
		self.name = name
		self.data = data

	@property
	def outpath(self):
		return "%s.svg" % self.name

	@property
	def style(self):
		return DefaultStyle(
			background='white',
			plot_background = 'white'
		)
	
	def line(self):
		pass

	def pie(self):
		chart = pygal.Pie(
			style = self.style,
			width = 550,
			height = 400,
			legend_at_bottom = True,
			legend_at_bottom_columns=6
		)
		for item, value in self.data:
			chart.add(item, value)
		chart.render_to_file(self.outpath)


class Statistics(object):
	meta_table = MetaTable()
	_bases = {}
	_total_sequences = None
	_total_bases = None
	unit = 'Mb'
	letter = 'ATGC'
	
	def __init__(self):
		self.table = FastaTable()
		
	def _calc_bases(self):
		'''
		calculate the number of bases A,T,G,C
		'''
		self._total_bases = 0
		for fasta_file in self.table.fetchAll():
			fastas = pyfaidx.Fasta(fasta_file.path, sequence_always_upper=True)
			for fasta in fastas:
				seq = fasta[:].seq
				self._total_bases += len(seq)

				for b in ['A','T','G','C']:
					self._bases[b] = self._bases.get(b, 0) + seq.count(b)

		for b in self._bases:
			self.meta_table.insert({'name': b, 'value': self._bases[b]})

		self.meta_table.insert({'name': 'total_bases', 'value': self._total_bases})

		return self._bases

	def _db_bases(self):
		'''
		get the number of bases A,T,G,C
		'''
		for data in self.meta_table.query("SELECT * FROM meta WHERE name IN ('A', 'T', 'G', 'C')"):
			self._bases[data.name] = data.value

		return self._bases

	@property
	def bases(self):
		if self._bases is None:
			self._db_bases() or self._calc_bases()

		return self._bases

	@property
	def seqcount(self):
		'''
		get the number of sequences in alayzed fastas
		'''
		if not self._total_sequences:
			seq_table = SequenceTable()
			self._total_sequences = seq_table.get("SELECT COUNT(1) FROM sequence LIMIT 1")

		return self._total_sequences

	@property
	def size(self):
		'''
		get sequence total length A+G+T+C+N
		'''
		if self._total_bases is None:
			self._total_bases = self.meta_table.getMeta('total_bases')

			if not self._total_bases:
				self._calc_bases()

		return self._total_bases

	@property
	def ns(self):
		'''
		get the number and percent of unknown bases
		@return (ns, percent)
		'''
		ns = self.size - self.validsize
		percent = round(ns/self.size*100, 2)
		return (ns, percent)

	@property
	def validsize(self):
		'''
		get sequence valid total length A+G+T+C
		'''
		return sum(self.bases.values())

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

class MicrosatelliteStatistics(Statistics):
	ssr_table = MicrosatelliteTable()
	_ssr_counts = None
	_ssr_length = None
	
	def __init__(self):
		super(MicrosatelliteStatistics, self).__init__()

	@property
	def counts(self):
		'''
		get total perfect SSR counts and frequency
		'''
		if self._ssr_counts is None:
			self._ssr_counts = self.ssr_table.recordCounts()
		
		return self._ssr_counts

	@property
	def lengths(self):
		if self._ssr_length is None:
			self._ssr_length = self.ssr_table.get("SELECT SUM(length) FROM ssr")
		
		return self._ssr_length

	@property
	def frequency(self):
		return self.ra(self.counts)

	@property
	def density(self):
		return self.rd(self.lengths)

	def getMotifLenStat(self):
		types = {1: 'Mono-', 2: 'Di-', 3: 'Tri-', 4: 'Tetra-', 5: 'Penta-', 6: 'Hexa-'}
		sql = "SELECT length(motif) AS type, SUM(length) AS length, COUNT(1) AS count FROM ssr GROUP BY type ORDER BY type"
		rows = [('Type', 'Counts', 'Length (bp)', 'Percent (%)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.ssr_table.query(sql):
			percent = round(row.count/self.counts*100, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((types[row.type], row.count, row.length, percent, frequency, density))
		return rows

	def getMotifTypeStat(self):
		sql = "SELECT standard, SUM(length) AS length, COUNT(1) AS count FROM ssr GROUP BY standard ORDER BY length(motif),standard"
		rows = [('Motif', 'Counts', 'Length (bp)', 'Percent (%)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.ssr_table.query(sql):
			percent = round(row.count/self.counts*100, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((row.standard, row.count, row.length, percent, frequency, density))
		return rows

	def getMotifRepeatStat(self):
		sql = "SELECT repeat, SUM(length) AS length, COUNT(1) AS count FROM ssr GROUP BY repeat ORDER BY repeat"
		rows = [('Repeat', 'Counts', 'Length (bp)', 'Percent (%)', 'Relative Abundance (loci/%s)' % self.unit, 'Relative Density (bp/%s)' % self.unit)]
		for row in self.ssr_table.query(sql):
			percent = round(row.count/self.counts*100, 2)
			frequency = self.ra(row.count)
			density = self.rd(row.length)
			rows.append((row.repeat, row.count, row.length, percent, frequency, density))
		return rows

class StatisticsReport:
	def __init__(self, unit='Mb', letter='ATGC'):
		self.data = Data()
		self.ssr_stat = MicrosatelliteStatistics()

	def generateDataTable(self):
		with open('ssr-statistics.json', 'wb') as fp:
			json.dump(self.data, fp)

	def readDataTable(self):
		with open('ssr-statistics.json', 'rb') as fh:
			data = json.load(fh)
		return Data(data)

	def sequenceSummary(self):
		self.data.seqcount = self.ssr_stat.seqcount
		self.data.size = self.ssr_stat.size
		self.data.validsize = self.ssr_stat.validsize
		self.data.ns = self.ssr_stat.ns
		self.data.gc = self.ssr_stat.gc
		self.data.unit = self.ssr_stat.unit

	def microsatelliteSummary(self):
		self.data.ssrcounts = self.ssr_stat.counts
		self.data.frequency = self.ssr_stat.frequency
		self.data.ssrlength = self.ssr_stat.lengths
		self.data.density = self.ssr_stat.density
		
		self.data.ssrtypes = self.ssr_stat.getMotifLenStat()
		chart_data = [row[0:2] for row in self.data.ssrtypes[1:]]
		Plots('niblet_ssr_len', chart_data).pie()

		self.data.ssrmotifs = self.ssr_stat.getMotifTypeStat()
		self.data.ssrrepeats = self.ssr_stat.getMotifRepeatStat()

