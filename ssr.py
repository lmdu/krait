#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import multiprocessing

import db
import motif
import zfasta
import tandem

class Task(multiprocessing.Process):
	deamon = True
	#def __init__(self):
		#super(Task, self).__init__()
		#self.deamon = True

	@property
	def db(self):
		if self._db is None:
			self._db = db.Database()
		return self._db

	def build_fasta_index(self, fasta_id, fasta_path):
		'''
		build index for fasta file and write fasta sequence to database
		@para fasta_id int, the fasta file id in database
		@para fasta_path str, the file path of fasta
		@return object, a Fasta object
		'''
		seqs = zfasta.Fasta(fasta_path)
		sql = "SELECT * FROM seq WHERE name='%s'" % seqs.keys()[0]
		if not self.db.get_one(sql):
			rows = [(None, name, fasta_id) for name in seqs.keys()]
			sql = "INSERT INTO seq VALUES (?,?,?)"
			self.db.get_cursor().executemany(sql, rows)

		return seqs

class SSRTask(Task):
	'''
	perfect microsatellite search thread
	'''
	def __init__(self, fastas, min_repeats, standard_level):
		super(SSRTask, self).__init__()
		self.fastas = fastas
		self.min_repeats = min_repeats
		self.motifs = motif.StandardMotif(standard_level)
		self.fasta_counts = len(self.fastas)
		self.progress = 0

	def run(self):
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			#insert ssr to database
			sql = "INSERT INTO ssr VALUES (?,?,?,?,?,?,?,?,?)"

			current_seqs = 0
			seq_counts = len(seqs)
			#start search perfect microsatellites
			for name, seq in seqs:
				current_seqs += 1
				seq_progress = current_seqs/seq_counts
				ssrs = tandem.search_ssr(seq, self.min_repeats)
				
				def values():
					for ssr in ssrs:
						row = [None, name, self.motifs.standard(ssr[0])]
						row.extend(ssr)
						yield row

				cursor = self.db.get_cursor()
				cursor.execute("BEGIN;")
				cursor.executemany(sql, values())
				cursor.execute("COMMIT;")
				self.progress = int(seq_progress*fasta_progress*100)

