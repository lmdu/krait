import csv
import time
import json
import numpy
import jinja2
import requests
import functools
import multiprocessing

from PySide2.QtCore import *

#import plot
import motif
from libs import *
from db import *
from utils import *
from config import *
from statistics import *

class Worker(QObject):
	update_progress = Signal(int)
	update_message = Signal(str)
	finished = Signal()
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
		seqs = fasta.GzipFasta(fasta_path)
		
		#get the first sequence name
		seqnames = seqs.keys()
		for first in seqnames:
			break

		sql = "SELECT 1 FROM seq WHERE name='{}' LIMIT 1".format(first)
		if not self.db.get_one(sql):
			rows = []
			for name in seqnames:
				row = (None, name, fasta_id, seqs.get_len(name), seqs.get_gc(name), seqs.get_ns(name))
				rows.append(row)
			self.db.insert("INSERT INTO seq VALUES (?,?,?,?,?,?)", rows)
		
		return seqs

	def emit_progress(self, percent):
		self.update_progress.emit(percent)

	def emit_message(self, msg):
		self.update_message.emit(msg)

	def emit_finish(self, msg):
		self.update_progress.emit(100)
		self.update_message.emit(msg)
		self.finished.emit()

	def process(self):
		pass

	def run(self):
		self.emit_progress(0)
		#try:
		self.process()
		#except Exception, e:
		#	self.emit_finish('Error: %s' % str(e))


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
			self.emit_message("Building fasta index for %s" % fasta_file)
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			total_bases = seqs.get_total_length()

			#insert ssr to database
			sql = "INSERT INTO ssr VALUES (?,?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			for name, seq in seqs:
				current_bases += seqs.get_len(name)
				seq_progress = current_bases/total_bases

				self.emit_message("Search perfect SSRs from %s" % name)
				ssrs = tandem.search_ssr(seq, self.min_repeats)

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
			self.emit_message("Building fasta index for %s" % fasta_file)

			seqs = self.build_fasta_index(fasta_id, fasta_file)
			total_bases = seqs.get_total_length()
			#insert ssr to database
			sql = "INSERT INTO issr VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			for name, seq in seqs:
				current_bases += seqs.get_len(name)
				seq_progress = current_bases/total_bases

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

	def process(self):
		self.db.set_option('cssr_start_time', int(time.time()))
		ssrs = self.db.query("SELECT * FROM ssr")
		total = self.db.get_one("SELECT COUNT(1) FROM ssr LIMIT 1")
		self.db.begin()
		self.emit_message("Concatenate compound SSRs...")
		cssrs = [ssrs.__next__()]
		prev_progress = 0

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

		self.db.commit()
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
		sql = "INSERT INTO cssr VALUES (?,?,?,?,?,?,?,?,?,?)"
		self.db.get_cursor().execute(sql,
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
			self.emit_message("Building fasta index for %s" % fasta_file)
			seqs = self.build_fasta_index(fasta_id, fasta_file)
			total_bases = seqs.get_total_length()
			#insert ssr to database
			sql = "INSERT INTO vntr VALUES (?,?,?,?,?,?,?,?)"

			current_bases = 0
			#start search perfect microsatellites
			for name, seq in seqs:
				current_bases += seqs.get_len(name)
				seq_progress = current_bases/total_bases

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
		selected = self.model.getSelectedRows()

		primerdesign.loadThermoParams(PRIMER3_CONFIG)
		primerdesign.setGlobals(self.primer3_settings, None, None)
		#total ssr counts in a table
		total_ssrs = self.db.get_one("SELECT COUNT(1) FROM %s" % table)
		total_select = len(selected)

		if total_ssrs == total_select:
			sql = "SELECT * FROM %s" % table
		else:
			sql = "SELECT * FROM %s WHERE id IN (%s) ORDER BY id" % (table, ",".join(map(str, selected)))

		def iter_ssrs():
			if total_ssrs == total_select:
				sql = "SELECT * FROM %s" % table
				for ssr in self.db.query(sql):
					yield ssr
			else:
				for sid in sorted(selected):
					sql = "SELECT * FROM %s WHERE id=%s" % (table, sid)
					yield self.db.get_row(sql)

		current = 0
		seqs = None
		succeeded = 0
		progress = 0
		prev_progress = 0

		insert_sql = "INSERT INTO primer VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
		
		self.db.begin()
		for item in iter_ssrs():
			if seqs is None or item.sequence not in seqs:
				sql = "SELECT f.path FROM fasta AS f,seq AS s WHERE f.id=s.fid AND s.name='{}' LIMIT 1".format(item.sequence)
				seqfile = self.db.get_one(sql)
				seqs = fasta.GzipFasta(seqfile)
			
			start = item.start - self.flank
			if start < 1:
				start = 1
			end = item.end + self.flank
			
			target = dict(
				SEQUENCE_ID = "%s-%s" % (table, item.id),
				SEQUENCE_TEMPLATE = seqs.get_seq_by_loci(item.sequence, start, end),
				SEQUENCE_TARGET = [item.start-start, item.length],
				SEQUENCE_INTERNAL_EXCLUDED_REGION = [item.start-start, item.length]
			)

			primerdesign.setSeqArgs(target)
			res = primerdesign.runDesign(False)
			current += 1

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

		self.db.commit()

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
		headers = self.model.columnNames()

		whole_counts = self.db.get_one("SELECT COUNT(*) FROM %s" % table_name)
		total_counts = whole_counts

		if self.selected == 'whole' or len(self.model.selected) == whole_counts:
			sql = "SELECT * FROM {}".format(table_name)
		else:
			ids = self.model.getSelectedRows()
			total_counts = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(ids))

		prev_progress = 0
		progress = 0
		current = 0

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

			for row in self.db.query(sql):
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
			ids = self.model.getSelectedRows()
			total_ssrs = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(ids))

		current = 0
		progress = 0
		prev_progress = 0
		seqs = None

		with open(self.outfile, 'wt') as fp:
			for item in self.db.query(sql):
				if seqs is None or item.sequence not in seqs:
					sql = "SELECT f.path FROM fasta AS f,seq AS s WHERE f.id=s.fid AND s.name='{}' LIMIT 1".format(item.sequence)
					seqfile = self.db.get_one(sql)
					seqs = fasta.GzipFasta(seqfile)
				
				start = item.start - self.flank
				if start < 1:
					start = 1
				end = item.end + self.flank
				ssr = seqs.get_seq_by_loci(item.sequence, start, end)
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

		table_name = self.model.tableName()
		headers = self.model.columnNames()

		whole_counts = self.db.get_one("SELECT COUNT(*) FROM %s" % table_name)
		total_counts = whole_counts

		if self.selected == 'whole' or len(self.model.selected) == whole_counts:
			sql = "SELECT * FROM {}".format(table_name)
		else:
			ids = self.model.getSelectedRows()
			total_counts = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(ids))

		prev_progress = 0
		progress = 0
		current = 0

		with open(self.outfile, 'w', newline='') as outfh:
			if self.outfile.endswith('.csv'):
				writer = csv.writer(outfh)
			else:
				writer = csv.writer(outfh, delimiter='\t')

			ssr_table = self.db.get_one("SELECT category FROM primer LIMIT 1")
			ssr_headers = self.db.get_fields(ssr_table)
			writer.writerow(ssr_headers + headers[3:])

			for row in self.db.query(sql):
				ssr = self.db.get_row("SELECT * FROM {} WHERE id={} LIMIT 1".format(row.category, row.target))

				print(row.category, row.target)

				ssr_with_primer = list(ssr.getValues())
				ssr_with_primer.extend(row.getValues()[3:])

				writer.writerow(ssr_with_primer)
				
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
		interval_forest = {}

		f = check_gene_annot_format(self.annot_file)

		if f == 'GFF':
			features = get_gff_coordinate(self.annot_file)
		else:
			features = get_gtf_coordinate(self.annot_file)

		locations = {}
		locid = 0
		prev_chrom = None
		starts = []
		ends = []
		indexes = []

		for feature in features:
			locid += 1
			locations[locid] = feature[3]

			if feature[0] != prev_chrom:
				if starts:
					starts = numpy.array(starts, dtype=numpy.long)
					ends = numpy.array(ends, dtype=numpy.long)
					indexes = numpy.array(indexes, dtype=numpy.long)
					interval_forest[prev_chrom] = ncls.NCLS(starts, ends, indexes)
				
				prev_chrom = feature[0]
				starts = []
				ends = []
				indexes = []

			starts.append(feature[1])
			ends.append(feature[2])
			indexes.append(locid)

		if starts:
			starts = numpy.array(starts, dtype=numpy.long)
			ends = numpy.array(ends, dtype=numpy.long)
			indexes = numpy.array(indexes, dtype=numpy.long)
			interval_forest[prev_chrom] = ncls.NCLS(starts, ends, indexes)

		if not interval_forest:
			self.emit_finish("No feature found in annotation file.")
			return

		#do mapping
		total = self.db.get_one("SELECT COUNT(1) FROM %s LIMIT 1" % self.table)
		current = 0
		progress = 0
		prev_progress = 0
		#categories = {'ssr': 1, 'cssr': 2, 'issr': 3, 'vntr': 4}
		#features = {'CDS': 1, 'UTR': 2, 'EXON': 3, 'INTRON': 4}
		#cat = categories[self.table]
		cat = self.table
		for ssr in self.db.get_cursor().execute("SELECT * FROM %s" % self.table):
			self.emit_message("Locating %ss in sequence %s" % (self.table.upper(), ssr.sequence))
			current += 1
			progress = int(current/total*100)
			if progress > prev_progress:
				self.emit_progress(progress)
				prev_progress = progress

			if ssr.sequence not in interval_forest:
				continue

			res = set(interval_forest[ssr.sequence].find_overlap(ssr.start, ssr.end))
			if not res:
				continue

			candidates = {'CDS', 'exon', 'intron', 'UTR', '5UTR', '3UTR'}
			feats = [locations[fid[2]] for fid in res if locations[fid[2]].feature in candidates]
			ssrfeats = set()
			for feat in sorted(feats, key=lambda x:x.feature):
				if feat.feature == 'exon':
					checks = [(i, feat.gene_id, feat.gene_name) in ssrfeats for i in ['CDS', '5UTR', '3UTR', 'UTR']]
					if any(checks):
						continue
				
				if (feat.feature, feat.gene_id, feat.gene_name)	not in ssrfeats:
					ssrfeats.add((feat.feature, feat.gene_id, feat.gene_name))

			for ssrfeat in ssrfeats:
				if cat == 'issr':
					repeats = round(ssr.length/ssr.type, 2)
				else:
					repeats = ssr.repeat
				record = [None, cat, ssr.id, ssr.motif, repeats, ssrfeat[0], ssrfeat[1], ssrfeat[2]]
				self.db.get_cursor().execute("INSERT INTO feature VALUES (?,?,?,?,?,?,?,?)", record)

		self.emit_finish("%s location completed." % self.table)

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
			ids = self.model.getSelectedRows()
			total_counts = len(ids)
			sql = "SELECT * FROM {} WHERE id IN ({})".format(table_name, ",".join(ids))

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
