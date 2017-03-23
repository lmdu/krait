#!/usr/bin/env python
# -*- coding: utf-8 -*-

from db import *
from fasta import *

class Detail(QObject):
	def __init__(self):
		self.flank = 100
		self.fasta = None
		self.db = Database()

	def getSequence(self, name, start, end):
		sequence = self.fasta.get_sequence_by_loci(name, start, end)
		left_start = start - self.flank
		if left_start < 1:
			left_start = 1
		left_flank = self.fasta.get_sequence_by_loci(name, left_start, end-1)
		right_flank = self.fasta.get_sequence_by_loci(name, end+1, end+self.flank)

		return sequence, left_flank, right_flank

	def formatBase(self, base):
		return '<span class="{0}">{0}</span>'.format(base)

	def formatNormal(self, bases):
		return "".join([self.formatBase(b) for b in bases])

	def formatTarget(self, target):
		return '<span class="target">%s</span>' % self.formatNormal(target)

	def generateHtml(self, body):
		html = """
		<html>
		<head>
		<style>
		.A,.T,.G,.C {
			display: inline-block;
			width: 24px;
			font-size: 16px;
			font-weight: bold;
			text-align: center;
		}
		.A {color: #5050ff;}
		.T {color: #e6e600;}
		.G {color: #00c000;}
		.C {color: #e00000;}
		.N {color: #7a7a7a;}

		.target span{border-bottom: 2px solid black;}
		</style>
		</head>
		<body>%s</body>
		</html>
		"""
		return html % body


class SSRDetail(Detail):
	def __init__(self, identify, flank=100):
		super(SSRDetail, self).__init__()
		self.flank = flank
		self.identify = identify

	def getSSR(self):
		sql = (
			"SELECT f.path,t.id,t.sequence,t.start,t.end,t.length FROM "
			"fasta AS f,seq AS s,ssr AS t WHERE f.id=s.fid AND "
			"t.sequence=s.name AND t.id=%s"
		)

		row = self.db.get_row(sql % self.identify)
		self.fasta = Fasta(row[0])
		seq, left, right = self.getSequence(row[2], row[3], row[4])
		content = "".join([self.formatNormal(left), self.formatTarget(seq), self.formatNormal(right)])

		return self.generateHtml(content)


class PrimerDetail(Detail):
	def __init__(self):
		super(PrimerDetail, self).__init__()
