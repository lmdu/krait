#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from PySide.QtCore import QThread, Signal
from ssr import *
from db import *
from utils import Data

class Worker(QThread):
	update_progress = Signal(int)
	update_message = Signal(str)
	
	def __init__(self, parent):
		super(Worker, self).__init__(parent)

class MicrosatelliteWorker(Worker):
	'''
	perfect microsatellite search thread
	'''
	def __init__(self, parent, fastas, rules):
		super(MicrosatelliteWorker, self).__init__(parent)
		self.fastas = fastas
		self.rules = rules
		self.fasta_counts = len(self.fastas)
		
		self.ssr_table = MicrosatelliteTable()
		self.seq_table = SequenceTable()

	def run(self):
		current_fastas = 0
		for fasta in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.update_message.emit("Building fasta index for %s." % fasta)
			detector = MicrosatelliteDetector(fasta.path, self.rules)

			#get all sequence names
			for name in detector.fastas.keys():
				self.seq_table.insert(Data(sid=None, name=name, fid=fasta.fid))

			#start search perfect SSRs
			self.update_message.emit("Searching perfect SSRs...")
			for ssr in detector:
				self.ssr_table.insert(ssr)
				self.update_progress.emit(round(fasta_progress*detector.progress*100))

			self.update_progress.emit(100)
			self.update_message.emit('Perfect SSRs search completed.')


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
			progress = round(int(cssr.cssrs.split(',')[-1])/total_ssrs*100)
			self.update_progress.emit(progress)

		self.update_progress.emit(100)
		self.update_message.emit("Compound SSRs search completed.")

class SatelliteWorker(Worker):
	def __init__(self, parent, fastas, motifs, repeats):
		super(SatelliteWorker, self).__init__(parent)
		self.fastas = fastas
		self.motifs = motifs
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
			detector = SatelliteDetector(fasta.path, self.motifs, self.repeats)

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

