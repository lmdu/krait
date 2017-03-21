#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import time
import json
import jinja2
import pyfaidx
from PySide.QtCore import *

import fasta
import motif
import tandem
import primerdesign
from ssr import *
from db import *
from fasta import *
from config import *
from statistics import StatisticsReport
from utils import Data

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

	def build_fasta_index(self, fasta_id, fasta_path):
		'''
		build index for fasta file and write fasta sequence to database
		@para fasta_id int, the fasta file id in database
		@para fasta_path str, the file path of fasta
		@return object, a Fasta object
		'''
		seqs = fasta.Fasta(fasta_path)
		sql = "SELECT * FROM seq WHERE name='%s'" % seqs.keys()[0]
		if not self.db.get_one(sql):
			rows = [(None, name, fasta_id) for name in seqs.keys()]
			sql = "INSERT INTO seq VALUES (?,?,?)"
			self.db.get_cursor().executemany(sql, rows)

		return seqs

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

	def process(self):
		start = time.time()
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.update_message.emit("Building fasta index for %s" % fasta_file)
			seqs = self.build_fasta_index(fasta_id, fasta_file)
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


class PrimerWorker(Worker):
	def __init__(self, table, ids, flank, primer3_settings):
		'''
		design primers for select row or all rows
		@para table str, the table names in database
		@para ids list, the selected row id
		@para flank int, the length of flanking sequence used to design primer
		@para primer3_settings str, the path of primer3 settings file
		'''
		super(PrimerWorker, self).__init__()
		self.table = table
		self.ids = ids
		self.flank = flank
		self.primer3_settings = primer3_settings

	def process(self):
		primerdesign.loadThermoParams(PRIMER3_CONFIG)
		primerdesign.setGlobals(self.primer3_settings, None, None)
		sql = (
			"SELECT f.path,t.id,t.sequence,t.start,t.end,t.length FROM "
			"fasta AS f,seq AS s,%s AS t WHERE f.id=s.fid AND "
			"t.sequence=s.name AND t.id IN (%s)"
		)
		sql = sql % (self.table, ",".join(map(str, self.ids)))

		total = len(self.ids)
		current = 0
		seqs = None

		insert_sql = "INSERT INTO primer VALUES (?,?,?,?,?,?,?,?,?)"
		cursor = self.db.get_cursor()
		cursor.execute("BEGIN TRANSACTION;")

		for item in self.db.get_cursor().execute(sql):
			if seqs is None or item[2] not in seqs:
				seqs = Fasta(item[0])
			
			start = item[3] - self.flank
			if start < 1:
				start = 1
			end = item[4] + self.flank
			
			target = {}
			target['SEQUENCE_ID'] = "%s-%s" % (self.table, item[1])
			target['SEQUENCE_TEMPLATE'] = seqs.get_sequence_by_loci(item[2], start, end)
			target['SEQUENCE_TARGET'] = [item[3]-start, item[5]]
			target['SEQUENCE_INTERNAL_EXCLUDED_REGION'] = target['SEQUENCE_TARGET']

			primerdesign.setSeqArgs(target)
			res = primerdesign.runDesign(False)
			current += 1

			primer_count = res['PRIMER_PAIR_NUM_RETURNED']
			for i in range(primer_count):
				primer = [None, target['SEQUENCE_ID'], i+1]
				primer.append(res['PRIMER_LEFT_%s_SEQUENCE' % i])
				primer.append(round(res['PRIMER_LEFT_%s_TM' % i], 2))
				primer.append(round(res['PRIMER_LEFT_%s_GC_PERCENT' % i], 2))
				primer.append(res['PRIMER_RIGHT_%s_SEQUENCE' % i])
				primer.append(round(res['PRIMER_RIGHT_%s_TM' % i], 2))
				primer.append(round(res['PRIMER_RIGHT_%s_GC_PERCENT' % i], 2))
				cursor.execute(insert_sql, primer)
		
			self.update_progress.emit(round(current/total*100))

		cursor.execute("COMMIT;")

		self.update_message.emit('Primer design completed')
		self.finished.emit()
