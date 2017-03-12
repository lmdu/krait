#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import time
import jinja2
import pyfaidx
from PySide.QtCore import *

import fasta
import motif
import tandem
from ssr import *
from db import *
from statistics import StatisticsReport
from utils import Data
from config import MAX_ROWS

class Worker(QObject):
	update_progress = Signal(int)
	update_message = Signal(str)
	finished = Signal()
	_db = None

	@property
	def db(self):
		if self._db is None:
			self._db = Database()
		return self._db

class SSRWorker(Worker):
	'''
	perfect microsatellite search thread
	'''
	def __init__(self, fastas, min_repeats, standard_level):
		super(SSRWorker, self).__init__()
		self.fastas = fastas
		self.min_repeats = min_repeats
		self.motifs = motif.StandardMotif(standard_level)
		self.fasta_counts = len(self.fastas)

	@Slot()
	def process(self):
		start = time.time()
		current_fastas = 0
		for _, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.update_message.emit("Building fasta index for %s" % fasta_file)
			seqs = fasta.Fasta(fasta_file)
			#insert ssr to database
			sql = "INSERT INTO ssr VALUES (?,?,?,?,?,?,?,?,?)"

			current_seqs = 0
			seq_counts = len(seqs)
			#start search perfect microsatellites
			for name, seq in seqs:
				self.update_message.emit("Search perfect SSRs from %s" % name)
				current_seqs += 1
				seq_progress = current_seqs/seq_counts
				ssrs = tandem.search_ssr(seq, self.min_repeats)
				
				def values():
					for ssr in ssrs:
						row = [None, name, self.motifs.standard(ssr[0])]
						row.extend(ssr)
						yield row

				self.update_message.emit("Insert SSRs to database")
				cursor = self.db.get_cursor()
				cursor.execute("BEGIN TRANSACTION;")
				cursor.executemany(sql, values())
				cursor.execute("COMMIT;")
				self.update_progress.emit(int(seq_progress*fasta_progress*100))
				del seq
		self.update_progress.emit(100)
		self.update_message.emit('Perfect SSRs search completed')
		self.finished.emit()
		print time.time() - start


class CompoundWorker(Worker):
	def __init__(self, parent, dmax):
		super(CompoundWorker, self).__init__(parent)
		self.dmax = dmax
		
		self.ssr_table = MicrosatelliteTable()
		self.cssr_table = CompoundTable()

	def run(self):
		total_ssrs = self.ssr_table.recordCounts()
		ssrs = self.ssr_table.fetchAll()
		self.update_message.emit("Search compound SSRs...")
		for cssr in CompoundDetector(ssrs, self.dmax):
			self.cssr_table.insert(cssr)
			progress = round(int(cssr.component.split(',')[-1])/total_ssrs*100)
			self.update_progress.emit(progress)

		self.update_progress.emit(100)
		self.update_message.emit("Compound SSRs search completed.")

class SatelliteWorker(Worker):
	def __init__(self, parent, fastas, min_motif, max_motif, repeats):
		super(SatelliteWorker, self).__init__(parent)
		self.fastas = fastas
		self.min_motif = min_motif
		self.max_motif = max_motif
		self.repeats = repeats
		self.fasta_counts = len(self.fastas)
		
		self.ssr_table = SatelliteTable()
		self.seq_table = SequenceTable()

	def run(self):
		current_fastas = 0
		for fasta in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.update_message.emit("Building fasta index for %s." % fasta)
			detector = SatelliteDetector(fasta.path, self.min_motif, self.max_motif, self.repeats)

			#get all sequence names
			for name in detector.fastas.keys():
				self.seq_table.insert(Data(sid=None, name=name, fid=fasta.fid))

			#start search perfect SSRs
			self.update_message.emit("Search satellites...")
			for ssr in detector:
				self.ssr_table.insert(ssr)
				self.update_progress.emit(round(fasta_progress*detector.progress*100))

			self.update_progress.emit(100)
			self.update_message.emit('Satellites search completed.')

class StatisticsWorker(Worker):
	def __init__(self, editor, parent):
		super(StatisticsWorker, self).__init__(parent)
		self.editor = editor

	def run(self):
		stat = StatisticsReport()
		stat.sequenceSummary()
		stat.microsatelliteSummary()
		stat.generateDataTable()
		data = stat.readDataTable()
		#change the rows of ssr motif statistics
		if data.ssrmotifs:
			rows = sorted(data.ssrmotifs[1:], key=lambda x: -x[1])
			data.ssrmotifs = [data.ssrmotifs[0]]
			data.ssrmotifs.extend(rows[0:MAX_ROWS])

		with open('report.html') as fh:
			content = jinja2.Template(fh.read()).render(**data)
		self.update_message.emit(content)

