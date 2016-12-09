#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from PySide.QtCore import QThread, Signal
from ssr import SSRDetector, CSSRDetector
from db import PerfectSSRTable, CompoundSSRTable, SequenceSSRTable
from utils import Data


class SSRSearchWorker(QThread):
	update_progress = Signal(int)
	update_message = Signal(str)

	def __init__(self, parent, fastas, rules):
		super(SSRSearchWorker, self).__init__(parent)
		self.fastas = fastas
		self.rules = rules
		self.ssr_table = PerfectSSRTable()
		self.seq_table = SequenceSSRTable()

	def run(self):
		current_fastas = 0
		for fasta in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/len(self.fastas)
			
			#use fasta and create fasta file index
			self.update_message.emit("Building fasta index for %s." % fasta)
			detector = SSRDetector(fasta.path, self.rules)

			#get all sequence names
			for name in detector.fastas.keys():
				self.seq_table.insert(Data(ID=None, name=name, fasta=fasta.ID))

			#start search perfect SSRs
			self.update_message.emit("Searching perfect SSRs...")
			for ssr in detector:
				self.ssr_table.insert(ssr)
				self.update_progress.emit(round(fasta_progress*detector.progress*100))

			self.update_progress.emit(100)
			self.update_message.emit('Perfect SSRs search completed.')


class CSSRSearchWorker(QThread):
	update_progress = Signal(int)
	update_message = Signal(str)
	
	def __init__(self, parent, dmax):
		super(CSSRSearchWorker, self).__init__(parent)
		self.dmax = dmax
		self.ssr_table = PerfectSSRTable()
		self.cssr_table = CompoundSSRTable()

	def run(self):
		total_ssrs = self.ssr_table.getTotalCounts()
		ssrs = self.ssr_table.fetchAll()
		for cssr in CSSRDetector(ssrs, self.dmax):
			self.cssr_table.insert(cssr)
			progress = round(int(cssr.cssrs.split(',')[-1])/total_ssrs*100)
			self.update_progress.emit(progress)

		self.update_progress.emit(100)
		self.update_message.emit("Compound SSRs search completed.")