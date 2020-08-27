import os
import csv
import time
import json
import numpy
import pyfastx
import requests
import traceback
import functools
import multiprocessing

from PySide2.QtCore import *
#from PySide.QtCore import *

import motif
from libs import *
from db import *
from gff import *
from utils import *
from config import *
from statistics import *
from config import ROOT_PATH

os.environ['PRIMER3HOME'] = ROOT_PATH
from primer3 import primerdesign

def build_full_index(fafile):
	_ = pyfastx.Fasta(fafile, full_index=True)

class Worker(QObject):
	update_progress = Signal(int)
	update_message = Signal(str)
	error_message = Signal(str)
	finished = Signal()
	failed = Signal()

	_db = None

	def __init__(self):
		super(Worker, self).__init__()

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
		@return Fasta object
		'''
		#seqs = fasta.GzipFasta(fasta_path)
		self.emit_message("Building fasta index for %s" % fasta_path)

		with multiprocessing.Pool() as pool:
			pool.apply_async(build_full_index, (fasta_path,))
			pool.close()
			pool.join()

		seqs = pyfastx.Fasta(fasta_path)
		
		#get sequence detail information
		#sql = "SELECT * FROM seq INNER JOIN fasta ON (seq.fid=fasta.id) WHERE fasta.path='{}' LIMIT 1".format(fasta_path)
		#if not self.db.get_one(sql):
		#	rows = []
		#	for seq in seqs:
		#		compos = seq.composition
		#		ns = sum(compos[b] for b in compos if b not in ['A', 'T', 'G', 'C'])
		#		row = (None, seq.name, fasta_id, len(seq), compos.get('G',0)+compos.get('C',0), ns)
		#		rows.append(row)
		#	self.db.insert("INSERT INTO seq VALUES (?,?,?,?,?,?)", rows)

		sql = "SELECT * FROM option WHERE name='gc_content'"
		if not self.db.get_one(sql):
			gc = seqs.gc_content
			compos = seqs.composition
			ns = sum(compos[b] for b in compos if b not in ['A', 'T', 'G', 'C'])
			self.db.insert("INSERT INTO option (name, value) VALUES (?,?)", [
				('total_base', str(seqs.size)),
				('total_seqs', str(len(seqs))),
				('gc_content', str(gc)),
				('unkown_base', str(ns))
			])

		self.total_bases = seqs.size
		seqs = pyfastx.Fasta(fasta_path, build_index=False)
		return seqs

	def emit_progress(self, percent):
		self.update_progress.emit(percent)

	def emit_message(self, msg):
		self.update_message.emit(msg)

	def emit_finish(self, msg):
		self.update_progress.emit(100)
		self.update_message.emit(msg)
		self.finished.emit()

	def emit_fail(self):
		self.update_progress.emit(100)
		self.failed.emit()

	def process(self):
		pass

	def run(self):
		self.emit_progress(0)
		try:
			self.process()
		except Exception:
			self.error_message.emit(traceback.format_exc())
		
		self.emit_fail()


class SSRWorker(Worker):
	"""
	perfect microsatellite search thread
	"""
	def __init__(self, fastas, min_repeats, standard_level):
		super(SSRWorker, self).__init__()
		self.fastas = fastas
		self.min_repeats = min_repeats
		self.motifs = motif.StandardMotif(standard_level)
		self.fasta_counts = len(self.fastas)
		self.progress = 0

		parameters = Data(
			mono = min_repeats[0],
			di = min_repeats[1],
			tri = min_repeats[2],
			tetra = min_repeats[3],
			penta = min_repeats[4],
			hexa = min_repeats[5],
			level = standard_level
		)
		self.db.set_option("ssr_parameters", json.dumps(parameters))

	def process(self):
		self.db.set_option('ssr_start_time', int(time.time()))
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			#total_bases = seqs.get_total_length()

			#insert ssr to database
			sql = "INSERT INTO ssr VALUES (?,?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			with multiprocessing.Pool() as pool:
				for name, seq in seqs:
				#for s in seqs:
				#	name, seq = s.name, s.seq

					current_bases += len(seq)

					seq_progress = current_bases/self.total_bases

					self.emit_message("Searching for perfect SSRs from %s" % name)
					#ssrs = tandem.search_ssr(seq, self.min_repeats)

					res = pool.apply_async(tandem.search_ssr, (seq, self.min_repeats))
					ssrs = res.get()

					def values():
						for ssr in ssrs:
							row = [None, name, self.motifs.standard(ssr[0])]
							row.extend(ssr)
							yield row

					self.db.insert(sql, values())

					self.emit_progress(int(seq_progress*fasta_progress*100))

		self.db.set_option('ssr_end_time', int(time.time()))
		self.emit_finish('Perfect SSRs search completed')


class ISSRWorker(Worker):
	'''
	perfect microsatellite search thread
	'''
	def __init__(self, fastas, seed_repeat, seed_length, max_edits, mis_penalty, gap_penalty, score, standard_level):
		super(ISSRWorker, self).__init__()
		self.fastas = fastas
		self.motifs = motif.StandardMotif(standard_level)
		self.fasta_counts = len(self.fastas)
		self.seed_repeat = seed_repeat
		self.seed_length = seed_length
		self.max_edits = max_edits
		self.mis_penalty = mis_penalty
		self.gap_penalty = gap_penalty
		self.score = score

		parameters = Data(
			seed_repeat = seed_repeat,
			seed_length = seed_length,
			max_edits = max_edits,
			mis_penalty = mis_penalty,
			gap_penalty = gap_penalty,
			min_score = score,
			level = standard_level
		)
		self.db.set_option('issr_parameters', json.dumps(parameters))

	def process(self):
		self.db.set_option('issr_start_time', int(time.time()))
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			#total_bases = seqs.get_total_length()
			#insert ssr to database
			sql = "INSERT INTO issr VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			for name, seq in seqs:
				current_bases += len(seq)
				seq_progress = current_bases/self.total_bases

				self.emit_message("Search imperfect SSRs from %s" % name)

				issrs = tandem.search_issr(seq, self.seed_repeat, self.seed_length, self.max_edits, self.mis_penalty, self.gap_penalty, self.score, 500)
				
				def values():
					for issr in issrs:
						row = [None, name, self.motifs.standard(issr[0])]
						row.extend(issr)
						yield row

				self.db.insert(sql, values())
				self.emit_progress(int(seq_progress*fasta_progress*100))
		self.db.set_option("issr_end_time", int(time.time()))
		self.emit_finish('Imperfect SSRs search completed')


class CSSRWorker(Worker):
	def __init__(self, dmax):
		super(CSSRWorker, self).__init__()
		self.dmax = dmax
		parameters = Data(dmax = dmax)
		self.db.set_option('cssr_parameters', json.dumps(parameters))
		self.sql = "INSERT INTO cssr VALUES (?,?,?,?,?,?,?,?,?,?)"

	def process(self):
		self.db.set_option('cssr_start_time', int(time.time()))
		ssrs = self.db.query("SELECT * FROM ssr")
		total = self.db.get_one("SELECT COUNT(1) FROM ssr LIMIT 1")
		self.emit_message("Concatenate compound SSRs...")
		cssrs = [next(ssrs)]
		prev_progress = 0
		#self.db.begin()
		for ssr in ssrs:
			d = ssr.start - cssrs[-1].end - 1
			if ssr.sequence == cssrs[-1].sequence and d <= self.dmax:
				cssrs.append(ssr)
			else:
				if len(cssrs) > 1:
					self.concatenate(cssrs)
					progress = int(cssrs[-1].id/total*100)
					if progress > prev_progress:
						self.emit_progress(progress)
						prev_progress = progress
				cssrs = [ssr]

		if len(cssrs) > 1:
			self.concatenate(cssrs)

		#self.db.commit()
		self.db.set_option('cssr_end_time', int(time.time()))
		self.emit_finish("Compound SSRs search completed")

	def concatenate(self, cssrs):
		seqname = cssrs[-1].sequence
		start = cssrs[0].start
		end = cssrs[-1].end
		complexity = len(cssrs)
		component = "%s-%s" % (cssrs[0].id, cssrs[-1].id)
		motif = "-".join([cssr.motif for cssr in cssrs])
		length = sum(cssr.length for cssr in cssrs)
		gap = sum(cssr.start-cssrs[idx].end-1 for idx, cssr in enumerate(cssrs[1:]))
		structure = "-".join(["(%s)%s" % (cssr.motif, cssr.repeat) for cssr in cssrs])
		self.db.get_cursor().execute(self.sql,
			(None, seqname, start, end, motif, complexity, length, gap, component, structure)
		)


class VNTRWorker(Worker):
	def __init__(self, fastas, min_motif, max_motif, repeats):
		super(VNTRWorker, self).__init__()
		self.fastas = fastas
		self.min_motif = min_motif
		self.max_motif = max_motif
		self.repeats = repeats
		self.fasta_counts = len(self.fastas)

		parameters = Data(
			min_motif = min_motif,
			max_motif = max_motif,
			min_repeat = repeats
		)

		self.db.set_option('vntr_parameters', json.dumps(parameters))

	def process(self):
		self.db.set_option('vntr_start_time', int(time.time()))
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			#total_bases = seqs.get_total_length()
			#insert ssr to database
			sql = "INSERT INTO vntr VALUES (?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			for name, seq in seqs:
				current_bases += len(seq)
				seq_progress = current_bases/self.total_bases

				self.emit_message("Search VNTRs from %s" % name)
				vntrs = tandem.search_vntr(seq, self.min_motif, self.max_motif, self.repeats)
				
				def values():
					for vntr in vntrs:
						row = [None, name]
						row.extend(vntr)
						yield row

				self.db.insert(sql, values())
				self.emit_progress(int(seq_progress*fasta_progress*100))
		self.db.set_option('vntr_end_time', int(time.time()))	
		self.emit_finish('VNTRs search completed')

class StatisWorker(Worker):
	def __init__(self, unit='Mb', letter='ATGC', dpi=300):
		super(StatisWorker, self).__init__()
		self.unit = unit
		self.letter = letter
		self.dpi = dpi

	def process(self):
		self.emit_message("Doing sequence statistics...")

		seq_statis = Statistics(self.unit, self.letter).results()

		self.db.set_option('seq_statis', json.dumps(seq_statis))

		if not self.db.is_empty('ssr'):
			self.emit_message("Doing perfect SSR statistics...")
			ssr_statis = SSRStatistics().results()
			self.db.set_option('ssr_statis', json.dumps(ssr_statis))
		else:
			self.db.set_option('ssr_statis', '[]')

		if not self.db.is_empty('issr'):
			self.emit_message("Doing imperfect SSR statistics...")
			issr_statis = ISSRStatistics().results()
			self.db.set_option('issr_statis', json.dumps(issr_statis))
		else:
			self.db.set_option('issr_statis', '[]')

		
		if not self.db.is_empty('cssr'):
			self.emit_message("Doing compound SSR statistics...")
			cssr_statis = CSSRStatistics().results()
			self.db.set_option('cssr_statis', json.dumps(cssr_statis))
		else:
			self.db.set_option('cssr_statis', '[]')


		if not self.db.is_empty('vntr'):
			self.emit_message("Doing VNTR statistics...")
			vntr_statis = VNTRStatistics().results()
			self.db.set_option('vntr_statis', json.dumps(vntr_statis))
		else:
			self.db.set_option('vntr_statis', '[]')

		self.emit_finish("Statistics was successfully completed")


class PrimerWorker(Worker):
	def __init__(self, model, flank, primer3_settings):
		'''
		design primers for select row or all rows
		@para model, table model
		@para flank int, the length of flanking sequence used to design primer
		@para primer3_settings str, the path of primer3 settings file
		'''
		super(PrimerWorker, self).__init__()
		self.model = model
		self.flank = flank
		self.primer3_settings = primer3_settings

	def process(self):
		self.emit_message("Designing primers...")
		table = self.model.tableName()
		selected = sorted(self.model.selected)

		primerdesign.loadThermoParams(PRIMER3_CONFIG)
		primerdesign.setGlobals(self.primer3_settings, None, None)
		#primer3.bindings.setP3Globals(self.primer3_settings)
		#total ssr counts in a table
		total_ssrs = self.db.get_one("SELECT COUNT(1) FROM %s" % table)
		total_select = len(selected)

		if total_ssrs == total_select:
			sql = "SELECT * FROM %s" % table
		else:
			sql = "SELECT * FROM %s WHERE id IN (%s) ORDER BY id" % (table, ",".join(map(str, selected)))

		#def iter_ssrs():
		#	if total_ssrs == total_select:
		#		sql = "SELECT * FROM %s" % table
		#		for ssr in self.db.query(sql):
		#			yield ssr
		#	else:
		#		for sid in sorted(selected):
		#			sql = "SELECT * FROM %s WHERE id=%s" % (table, sid)
		#			yield self.db.get_row(sql)

		current = 0
		current_seq = None
		current_name = None
		succeeded = 0
		progress = 0
		prev_progress = 0

		if not self.db.is_empty('primer'):
			self.db.get_cursor().execute("DELETE FROM primer")

		insert_sql = "INSERT INTO primer VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
		
		#self.db.begin()
		for item in self.db.query(sql):
			if item.sequence != current_name:
				sql = "SELECT f.path FROM fasta AS f,seq AS s WHERE f.id=s.fid AND s.name='{}' LIMIT 1".format(item.sequence)
				seqfile = self.db.get_one(sql)
				#seqs = fasta.GzipFasta(seqfile)
				seqs = pyfastx.Fasta(seqfile)
				current_seq = seqs[item.sequence].seq
				current_name = item.sequence
			
			start = item.start - self.flank

			if start < 1:
				start = 1
			end = item.end + self.flank
			
			target = dict(
				SEQUENCE_ID = "%s-%s" % (table, item.id),
				SEQUENCE_TEMPLATE = current_seq[start-1:end],
				SEQUENCE_TARGET = [item.start-start, item.length],
				SEQUENCE_INTERNAL_EXCLUDED_REGION = [item.start-start, item.length]
			)
			
			#primer3.bindings.setP3SeqArgs(target)
			#res = primer3.bindings.runP3Design(True)
			primerdesign.setSeqArgs(target)
			res = primerdesign.runDesign(False)
			
			current += 1
			
			if res is None:
				continue

			primer_count = res['PRIMER_PAIR_NUM_RETURNED']
			
			if primer_count:
				succeeded += 1
			
			for i in range(primer_count):
				primer = [None, table, item.id, i+1]
				primer.append(res['PRIMER_PAIR_%s_PRODUCT_SIZE' % i])
				primer.append(res['PRIMER_LEFT_%s_SEQUENCE' % i])
				primer.append(round(res['PRIMER_LEFT_%s_TM' % i], 2))
				primer.append(round(res['PRIMER_LEFT_%s_GC_PERCENT' % i], 2))
				primer.append(round(res['PRIMER_LEFT_%s_END_STABILITY' % i], 2))
				primer.append(res['PRIMER_RIGHT_%s_SEQUENCE' % i])
				primer.append(round(res['PRIMER_RIGHT_%s_TM' % i], 2))
				primer.append(round(res['PRIMER_RIGHT_%s_GC_PERCENT' % i], 2))
				primer.append(round(res['PRIMER_RIGHT_%s_END_STABILITY' % i], 2))

				self.db.get_cursor().execute(insert_sql, primer)

				meta = [self.db.get_last_insert_rowid()]
				meta.extend(res['PRIMER_LEFT_%s' % i])
				meta.extend(res['PRIMER_RIGHT_%s' % i])
				self.db.get_cursor().execute("INSERT INTO primer_meta VALUES (?,?,?,?,?)", meta)
			
			progress = int(current/total_select*100)
			if progress > prev_progress:
				self.emit_progress(progress)
				prev_progress = progress

		#self.db.commit()

		self.emit_finish('Primer design completed, %s succeed %s failed' % (succeeded, total_select-succeeded))


class ExportTableWorker(Worker):
	def __init__(self, selected, model, outfile):
		super(ExportTableWorker, self).__init__()
		self.selected = selected
		self.model = model
		self.outfile = outfile

	def process(self):
		#get selected ids from table model
		self.emit_message("Exporting to %s" % self.outfile)

		table_name = self.model.tableName()
		#headers = self.model.columnNames()

		whole_counts = self.db.get_one("SELECT COUNT(*) FROM %s LIMIT 1" % table_name)
		total_counts = whole_counts

		repeat_type = {'ssr':1, 'cssr':2, 'issr':3, 'vntr':4}.get(table_name)

		if self.selected == 'whole' or len(self.model.selected) == whole_counts:
			if self.db.is_empty('location'):
				sql = "SELECT * FROM {}".format(table_name)
			else:
				sql = (
					"select {0}.*,location.feature,gene.geneid,gene.genename,gene.biotype "
					"from {0} left join location on (location.target={0}.id) inner join gene "
					"on (gene.id=location.gid) WHERE location.reptype={1}"
				).format(table_name, repeat_type)
		else:
			ids = sorted(self.model.selected)
			total_counts = len(ids)

			if self.db.is_empty('location'):
				sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(map(str,ids)))
			else:
				sql = (
					"select {0}.*,location.feature,gene.geneid,gene.genename,gene.biotype "
					"from {0} left join location on (location.target={0}.id) inner join gene "
					"on (gene.id=location.gid) WHERE location.reptype={1} AND {0}.id IN ({2})"
				).format(table_name, repeat_type, ",".join(map(str,ids)))

		prev_progress = 0
		progress = 0
		current = 0

		cursor = self.db.get_cursor()

		if not self.db.is_empty('location'):
			def convert_feat(cursor, x):
				x = list(x)
				x[-4] = {1:'CDS', 2:'exon', 3:'3UTR', 4:'intron', 5:'5UTR'}.get(x[-4], '')
				return x

			cursor.setrowtrace(convert_feat)

		cursor.execute(sql)

		headers = [field[0] for field in cursor.getdescription()]

		with open(self.outfile, 'w', newline='') as outfh:
			if self.outfile.endswith('.csv'):
				writer = csv.writer(outfh)
			else:
				writer = csv.writer(outfh, delimiter='\t')

			if self.outfile.endswith('.gff'):
				writer.writerow(["##gff-version 3"])
				writer.writerow(["##generated by Krait {}".format(VERSION)])
				write_line = lambda x: writer.writerow(format_to_gff(table_name, x))
			else:
				writer.writerow(headers)
				write_line = lambda x: writer.writerow(x)

			#for row in self.db.query(sql):
			for row in cursor:
				write_line(row)
				current += 1
				process = int(current/total_counts*100)

				if process > prev_progress:
					self.emit_progress(process)
					prev_progress = process

		self.emit_finish("Successfully exported to %s" % self.outfile)


class ExportFastaWorker(Worker):
	def __init__(self, selected, model, flank, outfile):
		super(ExportFastaWorker, self).__init__()
		self.selected = selected
		self.model = model
		self.flank = flank
		self.outfile = outfile

	def process(self):
		self.emit_message("Exporting fasta sequence to %s" % self.outfile)
		
		table_name = self.model.tableName()
		
		whole_ssrs = self.db.get_one("SELECT COUNT(1) FROM %s" % table_name)
		total_ssrs = whole_ssrs

		if self.selected == 'whole' or len(self.model.selected) == whole_ssrs:
			sql = "SELECT * FROM {}".format(table_name)
		else:
			ids = sorted(self.model.selected)
			total_ssrs = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(map(str,ids)))

		current = 0
		progress = 0
		prev_progress = 0
		current_seq = None
		current_name = None

		with open(self.outfile, 'wt') as fp:
			for item in self.db.query(sql):
				if item.sequence != current_name:
					sql = "SELECT f.path FROM fasta AS f,seq AS s WHERE f.id=s.fid AND s.name='{}' LIMIT 1".format(item.sequence)
					seqfile = self.db.get_one(sql)
					seqs = pyfastx.Fasta(seqfile)
					current_seq = seqs[item.sequence].seq
					current_name = item.sequence
				
				start = item.start - self.flank
				if start < 1:
					start = 1
				end = item.end + self.flank
				#ssr = seqs.fetch(item.sequence, (start, end))
				ssr = current_seq[start-1:end]
				name = ">{}{} {}:{}-{}|motif:{}".format(table_name.upper(), item.id, item.sequence, item.start, item.end, item.motif)
				fp.write("{}\n{}".format(name, format_fasta_sequence(ssr, 70)))
				
				current += 1
				progress = int(current/total_ssrs*100)
				if progress > prev_progress:
					self.emit_progress(progress)
					prev_progress = progress

		self.emit_finish("Successfully exported to fasta %s" % self.outfile)

class ExportPrimerWorker(Worker):
	def __init__(self, selected, model, outfile):
		super(ExportPrimerWorker, self).__init__()
		self.selected = selected
		self.model = model
		self.outfile = outfile

	def process(self):
		#get selected ids from table model
		self.emit_message("Exporting to %s" % self.outfile)

		#table_name = self.model.tableName()
		#headers = self.model.columnNames()

		table_name = self.db.get_one("SELECT category FROM primer LIMIT 1")

		whole_counts = self.db.get_one("SELECT COUNT(*) FROM %s" % table_name)
		total_counts = whole_counts

		if self.selected == 'whole' or len(self.model.selected) == whole_counts:
			sql = (
				"SELECT primer.entry,primer.product,primer.forward,primer.tm1,"
				"primer.gc1,primer.stability1,primer.reverse,primer.tm2,primer.gc2,"
				"primer.stability2, {0}.* FROM primer left join {0} ON ({0}.id=primer.target) "
				"WHERE category='{0}'"
			).format(table_name)
		else:
			ids = sorted(self.model.selected)
			total_counts = len(ids)

			sql = (
				"SELECT primer.entry,primer.product,primer.forward,primer.tm1,"
				"primer.gc1,primer.stability1,primer.reverse,primer.tm2,primer.gc2,"
				"primer.stability2, {0}.* FROM primer left join {0} ON ({0}.id=primer.target) "
				"WHERE category='{0}' AND primer.id IN ({1})"
			).format(table_name, ",".join(map(str,ids)))

			#sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(map(str,ids)))

		prev_progress = 0
		progress = 0
		current = 0

		with open(self.outfile, 'w', newline='') as outfh:
			if self.outfile.endswith('.csv'):
				writer = csv.writer(outfh)
			else:
				writer = csv.writer(outfh, delimiter='\t')

			#ssr_table = self.db.get_one("SELECT category FROM primer LIMIT 1")
			#ssr_headers = self.db.get_fields(ssr_table)
			#writer.writerow(ssr_headers + headers[3:])
			
			cursor = self.db.get_cursor()
			cursor.execute(sql)

			headers = [field[0] for field in cursor.getdescription()]

			writer.writerow(headers[10:]+headers[:10])

			#for row in self.db.query(sql):
			#	ssr = self.db.get_row("SELECT * FROM {} WHERE id={} LIMIT 1".format(row.category, row.target))

			#	print(row.category, row.target)

			#	ssr_with_primer = list(ssr.getValues())
			#	ssr_with_primer.extend(row.getValues()[3:])

			#	writer.writerow(ssr_with_primer)

			for row in cursor:
				writer.writerow(row[10:]+row[:10])
				
				current += 1
				process = int(current/total_counts*100)

				if process > prev_progress:
					self.emit_progress(process)
					prev_progress = process

		self.emit_finish("Successfully exported to %s" % self.outfile)

class SaveProjectWorker(Worker):
	def __init__(self, dbfile):
		super(SaveProjectWorker, self).__init__()
		self.dbfile = dbfile

	def process(self):
		self.emit_message("Save project to %s" % self.dbfile)
		bak = self.db.save(self.dbfile)
		with bak as b:
			while not b.done:
				b.step(100)
				progress = int((b.pagecount-b.remaining)/b.pagecount*100)
				self.emit_progress(progress)
		self.emit_finish("Project has been successfully saved to %s" % self.dbfile)


class LocateWorker(Worker):
	"""
	Locate the SSRs in which region of genome
	@para table, the table name in database
	@para gene_annot, the genome annotation file, gff or gtf
	@para repeat_annot, the repeatmask output file contains TEs
	"""
	def __init__(self, table, annot_file=None):
		super(LocateWorker, self).__init__()
		self.table = table
		self.annot_file = annot_file

	def process(self):
		self.emit_message("Building interval tree for annotation...")

		_format = check_gene_annot_format(self.annot_file)

		if _format == 'GTF':
			mapper = GTFParser(self.annot_file)
		else:
			mapper = GFFParser(self.annot_file)

		if not self.db.is_empty('gene'):
			self.db.get_cursor().execute("DELETE FROM gene")

		self.db.get_cursor().executemany("INSERT INTO gene VALUES (?,?,?,?,?,?,?)", mapper.gene_info)

		#do mapping
		total = self.db.get_one("SELECT COUNT(*) FROM %s LIMIT 1" % self.table)
		current = 0
		progress = 0
		prev_progress = 0
		prev_seq = None

		repeat_types = {'ssr': 1, 'cssr': 2, 'issr': 3, 'vntr': 4}
		recs = []

		#remove previous results
		self.db.get_cursor().execute("DELETE FROM location WHERE reptype=%s" % repeat_types[self.table])

		for ssr in self.db.get_cursor().execute("SELECT * FROM %s" % self.table):
			if prev_seq != ssr.sequence:
				self.emit_message("Mapping %ss in sequence %s" % (self.table.upper(), ssr.sequence))
				prev_seq = ssr.sequence

			locs = mapper.mapping(ssr.sequence, ssr.start, ssr.end)
			if locs is None:
				continue

			rec = [None, repeat_types[self.table], ssr.id, locs[0], locs[1]]
			recs.append(rec)

			current += 1
			progress = int(current/total*100)
			if progress > prev_progress:
				self.emit_progress(progress)
				prev_progress = progress

		self.db.get_cursor().executemany("INSERT INTO location VALUES (?,?,?,?,?)", recs)
		self.emit_finish("%s mapping completed." % self.table)

class ExportFeatureWorker(Worker):
	def __init__(self, selected, model, outfile):
		super(ExportFeatureWorker, self).__init__()
		self.selected = selected
		self.model = model
		self.outfile = outfile

	def process(self):
		#get selected ids from table model
		self.emit_message("Exporting to %s" % self.outfile)

		table_name = self.model.tableName()
		headers = self.model.columnNames()

		whole_counts = self.db.get_one("SELECT COUNT(*) FROM %s" % table_name)
		total_counts = whole_counts

		if self.selected == 'whole' or len(self.model.selected) == whole_counts:
			sql = "SELECT * FROM {}".format(table_name)
		else:
			ids = sorted(self.model.selected)
			total_counts = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(map(str,ids)))

		prev_progress = 0
		progress = 0
		current = 0

		ssr_table = self.db.get_one("SELECT category FROM feature LIMIT 1")
		ssr_headers = self.db.get_fields(ssr_table)

		with open(self.outfile, 'w', newline='') as outfh:
			if self.outfile.endswith('.csv'):
				writer = csv.writer(outfh)
			else:
				writer = csv.writer(outfh, delimiter='\t')

			writer.writerow(ssr_headers + ["GI;GN;LOC"])

			prev_ssr_id = None
			ssr_group = []
			for row in self.db.query(sql):
				if row.target != prev_ssr_id:
					if ssr_group:
						locations = ['{};{};{}'.format(marker.geneid, marker.genename, marker.location) for marker in ssr_group]
						ssr = self.db.get_row("SELECT * FROM {} WHERE id={} LIMIT 1".format(ssr_table, prev_ssr_id))
						ssr_with_primer = list(ssr.getValues())
						ssr_with_primer.extend(locations)
						writer.writerow(ssr_with_primer)
					
					ssr_group = [row]
					prev_ssr_id = row.target
				
				else:
					ssr_group.append(row)

				current += 1
				progress = int(current/total_counts*100)

				if progress > prev_progress:
					self.emit_progress(progress)
					prev_progress = progress

			if ssr_group:
				locations = ['{};{};{}'.format(ssr.geneid, ssr.genename, ssr.location) for ssr in ssr_group]
				ssr = self.db.get_row("SELECT * FROM {} WHERE id={} LIMIT 1".format(ssr_table, prev_ssr_id))
				ssr_with_primer = list(ssr.getValues())
				ssr_with_primer.extend(locations)
				writer.writerow(ssr_with_primer)

		self.emit_finish("Successfully exported to %s" % self.outfile)

class EutilWorker(Worker):
	def __init__(self, acc, outfile, bank='nucleotide'):
		super(EutilWorker, self).__init__()
		self.acc = acc
		self.outfile = outfile
		self.bank = bank
		self.base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=%s&rettype=fasta&id=%s'

	def process(self):
		url = self.base % (self.bank, self.acc)
		self.emit_message("Downloading %s fasta sequence from NCBI..." % self.acc)
		r = requests.get(url, timeout=10, stream=True)
		if r.status_code == requests.codes.ok:
			self.total = 0
			self.start = time.time()
			time.sleep(0.1)
			with open(self.outfile, "wb") as fh:
				for chunk in r.iter_content(chunk_size=1024):
					self.total += len(chunk)
					fh.write(chunk)
					self.emit_message(self.progressing())

			self.emit_finish("Download %s fasta completed" % self.acc)	
		else:
			self.emit_finish("%s error: %s" % (r.status_code, r.reason))

	def progressing(self):
		total = human_size(self.total)
		timer = time.time() - self.start
		speed = human_size(self.total/timer)
		return "Downloaded %s, Speed %s/s" % (total, speed)
