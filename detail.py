#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db import *
from fasta import *
from utils import *

class Detail(object):
	def __init__(self):
		self.fasta = None
		self.db = Database()

	def getSequence(self, name, start, end):
		'''
		get the sequence of tandem repeat (TR) with left and right
		flanking sequences
		@para name str, the sequence name where TR located
		@para start int, the 1-based start location of TR
		@para end int, the 1-based end location of TR
		@return tuple, (left flank, self, right flank)
		'''
		sequence = self.fasta.get_sequence_by_loci(name, start, end)
		left_start = start - self.flank
		if left_start < 1:
			left_start = 1
		left_flank = self.fasta.get_sequence_by_loci(name, left_start, start-1)
		right_flank = self.fasta.get_sequence_by_loci(name, end+1, end+self.flank)

		return sequence, left_flank, right_flank

	def formatTarget(self, bases):
		return "".join('<span class="{0}">{0}</span>'.format(b) for b in bases)

	def formatPrimer(self, flank, start, length):
		res = []
		for i, b in enumerate(flank):
			if start < i < start+length:
				res.append('<span class="{0}">{0}</span>'.format(b))
			else:
				res.append(b)

		return "".join(res)

class SequenceDetail(Detail):
	def __init__(self, table, identify, flank=100):
		super(SequenceDetail, self).__init__()
		self.table = table
		self.id = identify
		self.flank = flank

	def generateHtml(self):
		sql = (
			"SELECT f.path FROM fasta AS f, seq AS s, {0} AS t "
			"WHERE f.id=s.fid AND t.sequence=s.name AND t.id={1}"
		)
		fasta_file = self.db.get_one(sql.format(self.table, self.id))
		self.fasta = Fasta(fasta_file)

		sql = "SELECT * FROM %s WHERE id=%s" % (self.table, self.id)
		ssr = self.db.get_row(sql)
		seq, left, right = self.getSequence(ssr.sequence, ssr.start, ssr.end)

		tandem = "%s%s%s" % (left, self.formatTarget(seq), right)

		return template_render("sequence.html", tandem=tandem, ssr=ssr)


class PrimerDetail(Detail):
	def __init__(self, table, identify, flank):
		super(PrimerDetail, self).__init__()
		self.table = table
		self.id = identify
		self.flank = flank

	def generateHtml(self):
		sql = "SELECT * FROM primer,primer_meta WHERE id=pid AND id=%s" % self.id
		primer = self.db.get_row(sql)

		table, tid = primer.target.split('-')

		sql = (
			"SELECT f.path FROM fasta AS f, seq AS s, {0} AS t "
			"WHERE f.id=s.fid AND t.sequence=s.name AND t.id={1}"
		)

		fasta_file = self.db.get_one(sql.format(table, tid))
		self.fasta = Fasta(fasta_file)

		sql = "SELECT * FROM %s WHERE id=%s" % (table, tid)
		ssr = self.db.get_row(sql)
		seq, left, right = self.getSequence(ssr.sequence, ssr.start, ssr.end)

		tandem = "%s%s%s" % (
			self.formatPrimer(left, primer.start1, primer.length1),
			self.formatTarget(seq),
			self.formatPrimer(right, primer.start2, primer.length2)
		)

		return template_render("sequence.html", tandem=tandem, ssr=ssr)
