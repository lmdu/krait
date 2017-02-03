#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from PySide.QtSql import *
import pyfaidx
from db import *

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
		return round(gc/self.validsize, 2)

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
	html = '''<html><head><style>
	*{font-size: 13px;}
	td,th{padding: 5px 10px;}
	</style></head><body>%s</body></html>'''
	
	def __init__(self, unit='Mb', letter='ATGC'):
		self.contents = []
		self.stat = Statistics()
		self.unit = unit
		self.stat.unit = unit
		self.stat.letter = letter
		
		self.ssr_stat = MicrosatelliteStatistics()
	
	def table(self, rows):
		body = ['<table border="0" bgcolor="gray" cellpadding="0" cellspacing="1" align="center" width="80%">']
		body.append('<tr bgcolor="#ffffff">%s</tr>' % "".join(['<th>%s</th>' % title for title in rows[0]]))
		row_format = '<td>%s</td>' * len(rows[0])
		for row in rows[1:]:
			body.append('<tr bgcolor="#ffffff">%s</tr>' % row_format % row)
		body.append('</table>')
		return "".join(body)

	def add(self, content, alignment=None):
		if alignment:
			self.contents.append('<p style="text-align:%s">%s</p>' % (alignment, content))
		else:
			self.contents.append('<p>%s</p>' % content)

	def generateReport(self):
		return self.html % "".join(self.contents)

	def sequenceSummary(self):
		self.add("Total number of sequences: <b>%s</b>" % self.stat.seqcount)
		self.add("Total length of sequences (A+T+C+G+N): <b>%s</b> bp" % self.stat.size)
		self.add("Total valid length of sequences (A+T+C+G): <b>%s</b> bp" % self.stat.validsize)
		self.add("Unkown bases (Ns): <b>%s<b> (<b>%s%%</b>)" % self.stat.ns)
		self.add("GC content (G+C)/(A+T+C+G) not include Ns: <b>%s%%</b>" % self.stat.gc)

	def microsatelliteSummary(self):
		_ = "Total number of SSRs: <b>%s</b>, relative abundance: <b>%s</b> loci/%s"
		self.add(_ % (self.ssr_stat.counts, self.ssr_stat.frequency, self.unit))
		_ = "Total length of SSRs: <b>%s</b> bp, relative density: <b>%s</b> bp/%s"
		self.add(_ % (self.ssr_stat.lengths, self.ssr_stat.density, self.unit))
		self.add("The motif length summary:")
		self.add(self.table(self.ssr_stat.getMotifLenStat()), alignment='center')
		self.add("The motif type summary:")
		self.add(self.table(self.ssr_stat.getMotifTypeStat()), alignment='center')
		self.add("The motif repeat summary:")
		self.add(self.table(self.ssr_stat.getMotifRepeatStat()), alignment='center')




