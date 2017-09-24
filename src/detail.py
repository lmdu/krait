#!/usr/bin/env python
# -*- coding: utf-8 -*-
from libs import *
from db import *
from utils import *

class Detail(object):
	def __init__(self):
		self.fasta = None
		self.db = Database()

	def getSequence(self, chrom, start, end):
		'''
		get the sequence of tandem repeat (TR) with left and right
		flanking sequences
		@para chrom str, the sequence name where TR located
		@para start int, the 1-based start location of TR
		@para end int, the 1-based end location of TR
		@return tuple, (left flank, self, right flank)
		'''
		sequence = self.fasta[chrom][start-1:end]
		left_start = start - self.flank
		
		if left_start < 1:
			left_start = 1

		lflank = self.fasta[chrom][left_start:start-1]
		rflank = self.fasta[chrom][end:end+self.flank]

		return sequence, lflank, rflank

	def formatBase(self, b, target=False):
		if target:
			return '<div class="base {0}">{0}</div>'.format(b)
		else:
			return '<div class="base">%s</div>' % b

	def formatTarget(self, bases):
		return "".join(self.formatBase(b, True) for b in bases)

	def formatFlank(self, bases):
		return "".join(self.formatBase(b) for b in bases)

	def formatPrimer(self, flank, start, length):
		res = []
		for i, b in enumerate(flank):
			if start <= i+1 <= start+length-1:
				res.append(self.formatBase(b, True))
			else:
				res.append(self.formatBase(b))
		return "".join(res)

class SequenceDetail(Detail):
	def __init__(self, table, _id, flank=100):
		super(SequenceDetail, self).__init__()
		self.table = table
		self.id = _id
		self.flank = flank

	def generateHtml(self):
		sql = (
			"SELECT f.path FROM fasta AS f, seq AS s, {0} AS t "
			"WHERE f.id=s.fid AND t.sequence=s.name AND t.id={1}"
		)
		fasta_file = self.db.get_one(sql.format(self.table, self.id))
		self.fasta = fasta.GzipFasta(fasta_file)

		sql = "SELECT * FROM %s WHERE id=%s" % (self.table, self.id)
		ssr = self.db.get_row(sql)
		seq, left, right = self.getSequence(ssr.sequence, ssr.start, ssr.end)

		tandem = "%s%s%s" % (self.formatFlank(left), self.formatTarget(seq), self.formatFlank(right))

		return template_render("sequence.html", tandem=tandem, ssr=ssr, table=self.table)


class PrimerDetail(Detail):
	def __init__(self, table, _id, flank):
		super(PrimerDetail, self).__init__()
		self.table = table
		self.id = _id
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
		self.fasta = fasta.GzipFasta(fasta_file)

		sql = "SELECT * FROM %s WHERE id=%s" % (table, tid)
		ssr = self.db.get_row(sql)
		seq, left, right = self.getSequence(ssr.sequence, ssr.start, ssr.end)

		tandem = "%s%s%s" % (
			self.formatPrimer(left, primer.start1, primer.length1),
			self.formatTarget(seq),
			self.formatPrimer(right, primer.start2-primer.length2-len(seq)-len(left)+1, primer.length2)
		)

		return template_render("sequence.html", tandem=tandem, ssr=ssr, table=self.table)
