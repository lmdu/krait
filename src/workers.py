#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import time
import json
import jinja2
import requests
from PySide.QtCore import *

import plot
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
		self.emit_progress(0)

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
		@return tuple, (Fasta object, sequence counts)
		'''
		seqs = fasta.GzipFasta(fasta_path)
		sql = "SELECT 1 FROM seq WHERE name='%s' LIMIT 1" % seqs.keys[0]
		if not self.db.get_one(sql):
			rows = [(None, name, fasta_id) for name in seqs.keys]
			sql = "INSERT INTO seq VALUES (?,?,?)"
			self.db.insert(sql, rows)

		return seqs, len(seqs.keys)

	def emit_progress(self, percent):
		self.update_progress.emit(percent)
		time.sleep(0)

	def emit_message(self, msg):
		self.update_message.emit(msg)
		time.sleep(0)

	def emit_finish(self, msg):
		self.update_progress.emit(100)
		self.update_message.emit(msg)
		self.finished.emit()


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

	def process(self):
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.emit_message("Building fasta index for %s" % fasta_file)

			seqs, seq_counts = self.build_fasta_index(fasta_id, fasta_file)
			#insert ssr to database
			sql = "INSERT INTO ssr VALUES (?,?,?,?,?,?,?,?,?)"

			current_seqs = 0
			#start search perfect microsatellites
			with seqs:
				for name, seq in seqs:
					current_seqs += 1
					seq_progress = current_seqs/seq_counts

					self.emit_message("Search perfect SSRs from %s" % name)
					ssrs = tandem.search_ssr(seq, self.min_repeats)
					
					def values():
						for ssr in ssrs:
							row = [None, name, self.motifs.standard(ssr[0])]
							row.extend(ssr)
							yield row

					self.db.insert(sql, values())
					self.emit_progress(int(seq_progress*fasta_progress*100))

		self.emit_finish('Perfect SSRs search completed')


class ISSRWorker(Worker):
	'''
	perfect microsatellite search thread
	'''
	def __init__(self, fastas, seed_repeat, seed_length, max_eidts, mis_penalty, gap_penalty, score, standard_level):
		super(ISSRWorker, self).__init__()
		self.fastas = fastas
		self.motifs = motif.StandardMotif(standard_level)
		self.fasta_counts = len(self.fastas)
		self.seed_repeat = seed_repeat
		self.seed_length = seed_length
		self.max_eidts = max_eidts
		self.mis_penalty = mis_penalty
		self.gap_penalty = gap_penalty
		self.score = score

	def process(self):
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.emit_message("Building fasta index for %s" % fasta_file)

			seqs, seq_counts = self.build_fasta_index(fasta_id, fasta_file)
			#insert ssr to database
			sql = "INSERT INTO issr VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"

			current_seqs = 0
			#start search perfect microsatellites
			with seqs:
				for name, seq in seqs:
					current_seqs += 1
					seq_progress = current_seqs/seq_counts

					self.emit_message("Search imperfect SSRs from %s" % name)

					issrs = tandem.search_issr(seq, self.seed_repeat, self.seed_length, self.max_eidts, self.mis_penalty, self.gap_penalty, self.score, 2000)
					
					def values():
						for issr in issrs:
							row = [None, name, self.motifs.standard(issr[0])]
							row.extend(issr)
							yield row

					self.db.insert(sql, values())
					self.emit_progress(int(seq_progress*fasta_progress*100))
	
		self.emit_finish('Imperfect SSRs search completed')


class CSSRWorker(Worker):
	def __init__(self, dmax):
		super(CSSRWorker, self).__init__()
		self.dmax = dmax

	def process(self):
		ssrs = self.db.query("SELECT * FROM ssr")
		total = self.db.get_one("SELECT COUNT(1) FROM ssr LIMIT 1")
		self.db.begin()
		self.emit_message("Concatenate compound SSRs")
		cssrs = [ssrs.next()]
		for ssr in ssrs:
			d = ssr.start - cssrs[-1].end - 1
			if ssr.sequence == cssrs[-1].sequence and d <= self.dmax:
				cssrs.append(ssr)
			else:
				if len(cssrs) > 1:
					self.concatenate(cssrs)
					progress = int(cssrs[-1].id/total*100)
					#self.emit_progress(progress)
				cssrs = [ssr]

		if len(cssrs) > 1:
			self.concatenate(cssrs)

		self.db.commit()

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

	def process(self):
		current_fastas = 0
		for fasta_id, fasta_file in self.fastas:
			current_fastas += 1
			fasta_progress = current_fastas/self.fasta_counts
			
			#use fasta and create fasta file index
			self.emit_message("Building fasta index for %s" % fasta_file)
			seqs, seq_counts = self.build_fasta_index(fasta_id, fasta_file)
			#insert ssr to database
			sql = "INSERT INTO vntr VALUES (?,?,?,?,?,?,?,?)"

			current_seqs = 0
			#start search perfect microsatellites
			with seqs:
				for name, seq in seqs:
					current_seqs += 1
					seq_progress = current_seqs/seq_counts

					self.emit_message("Search VNTRs from %s" % name)
					vntrs = tandem.search_vntr(seq, self.min_motif, self.max_motif, self.repeats)
					
					def values():
						for vntr in vntrs:
							row = [None, name]
							row.extend(vntr)
							yield row

					self.db.insert(sql, values())
					self.emit_progress(int(seq_progress*fasta_progress*100)) 
				
		self.emit_finish('VNTRs search completed')

class StatisWorker(Worker):
	def __init__(self, unit='Mb', letter='ATGC', dpi=300):
		super(StatisWorker, self).__init__()
		self.unit = unit
		self.letter = letter
		self.dpi = dpi

	def process(self):
		self.emit_message("Doing sequence statistics...")

		try:
			seq_statis = Statistics(self.unit, self.letter).results()
		except Exception, e:
			self.emit_finish(str(e))
			return

		self.db.set_option('seq_statis', json.dumps(seq_statis))

		if not self.db.is_empty('ssr'):
			self.emit_message("Doing perfect SSR statistics...")
			ssr_statis = SSRStatistics().results()
			self.db.set_option('ssr_statis', json.dumps(ssr_statis))

			#generate ssr type distribution pie plot
			x = [row[1] for row in ssr_statis.type[1:]]
			l = [row[0] for row in ssr_statis.type[1:]]
			plot.pie(x, l, "ssr_type", self.dpi)

			#generate most abundant ssr motif distribution bar plot
			motifs = [[], [], [], [], [], []]
			for row in ssr_statis.category[1:]:
				motifs[len(row[0])-1].append((row[0], row[1]))

			x = []
			l1 = []
			for m in motifs:
				m = sorted(m, key=lambda x: (x[0], -x[1]))
				for a, b in m[:10]:
					x.append(b)
					l1.append(a)

			plot.bar(l1, x, "SSR motif category", "SSR counts", "ssr_motif", self.dpi)

			#generate ssr repeat distribution box plot
			x = ssr_statis.repeat
			plot.box(x, l, "SSR repeats", "ssr_repeat", self.dpi)

			#generate ssr length distribution box plot
			x = ssr_statis.ssrlen
			plot.box(x, l, "SSR length (bp)", "ssr_length", self.dpi)


			#generate ssr distribution in diff regions pie plot
			if ssr_statis.region:
				x = [row[1] for row in ssr_statis.region]
				l = [row[0] for row in ssr_statis.region]
				plot.pie(x, l, "ssr_region", self.dpi)

		else:
			self.db.set_option('ssr_statis', '[]')


		if not self.db.is_empty('issr'):
			self.emit_message("Doing imperfect SSR statistics...")
			issr_statis = ISSRStatistics().results()
			self.db.set_option('issr_statis', json.dumps(issr_statis))

			#generate issr type distribution pie plot
			x = [row[1] for row in issr_statis.type[1:]]
			l = [row[0] for row in issr_statis.type[1:]]
			plot.pie(x, l, "issr_type", self.dpi)

			#generate ssr repeat distribution box plot
			x = issr_statis.score
			plot.box(x, l, "iSSR score", "issr_score", self.dpi)

			#generate ssr length distribution box plot
			x = issr_statis.issrlen
			plot.box(x, l, "iSSR length (bp)", "issr_length", self.dpi)

			#generate ssr distribution in diff regions pie plot
			if issr_statis.region:
				x = [row[1] for row in issr_statis.region]
				l = [row[0] for row in issr_statis.region]
				plot.pie(x, l, "issr_region", self.dpi)

		else:
			self.db.set_option('issr_statis', '[]')

		
		if not self.db.is_empty('cssr'):
			self.emit_message("Doing compound SSR statistics...")
			cssr_statis = CSSRStatistics().results()
			self.db.set_option('cssr_statis', json.dumps(cssr_statis))

			#generate cssr complexity distribution
			x = [row[0] for row in cssr_statis.complexity[1:]]
			y = [row[1] for row in cssr_statis.complexity[1:]]
			plot.line(x, y, 'cSSR complexity', 'cSSR Counts', 'cssr_complexity')

			#genrate cssr length distribution
			x = [row[0] for row in cssr_statis.cssrlen[1:]]
			y = [row[1] for row in cssr_statis.cssrlen[1:]]
			plot.line(x, y, 'cSSR length (bp)', 'cSSR Counts', 'cssr_length')

			#genrate cssr gap distribution
			x = [row[0] for row in cssr_statis.gap[1:]]
			y = [row[1] for row in cssr_statis.gap[1:]]
			plot.line(x, y, 'Gap length (bp)', 'cSSR Counts', 'cssr_gap')

			#generate ssr distribution in diff regions pie plot
			if cssr_statis.region:
				x = [row[1] for row in cssr_statis.region]
				l = [row[0] for row in cssr_statis.region]
				plot.pie(x, l, "cssr_region", self.dpi)

		else:
			self.db.set_option('cssr_statis', '[]')


		if not self.db.is_empty('vntr'):
			self.emit_message("Doing VNTR statistics...")
			vntr_statis = VNTRStatistics().results()
			self.db.set_option('vntr_statis', json.dumps(vntr_statis))

			#generate vntr type distribution
			x = [row[0] for row in vntr_statis.type]
			y = [row[1] for row in vntr_statis.type]
			plot.line(x, y, 'VNTR motif length (bp)', 'VNTR Counts', 'vntr_type', self.dpi)

			#genrate vntr length distribution
			x = [row[0] for row in vntr_statis.vntrlen]
			y = [row[1] for row in vntr_statis.vntrlen]
			plot.line(x, y, 'VNTR length (bp)', 'cSSR Counts', 'vntr_length', self.dpi)

			#genrate vntr repeat distribution
			x = [row[0] for row in vntr_statis.repeat]
			y = [row[1] for row in vntr_statis.repeat]
			plot.line(x, y, 'VNTR repeats', 'cSSR Counts', 'vntr_repeat', self.dpi)

			#generate ssr distribution in diff regions pie plot
			if vntr_statis.region:
				x = [row[1] for row in vntr_statis.region]
				l = [row[0] for row in vntr_statis.region]
				plot.pie(x, l, "vntr_region", self.dpi)

		else:
			self.db.set_option('vntr_statis', '[]')

		self.emit_finish("Statistics was completed")


class PrimerWorker(Worker):
	def __init__(self, table, ids, flank, primer3_settings):
		'''
		design primers for select row or all rows
		@para table str, the table names in database
		@para ids dict, the selected row id
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
		#total ssr counts in a table
		total_ssrs = self.db.get_one("SELECT COUNT(1) FROM %s" % self.table)
		total_ids = len(self.ids)

		if total_ssrs == total_ids:
			sql = (
				"SELECT f.path,t.id,t.sequence,t.start,t.end,t.length FROM "
				"%s AS t,fasta AS f,seq AS s WHERE f.id=s.fid AND "
				"t.sequence=s.name"
			)
			sql = sql % self.table
		else:
			sql = (
				"SELECT f.path,t.id,t.sequence,t.start,t.end,t.length FROM "
				"%s AS t,fasta AS f,seq AS s WHERE f.id=s.fid AND "
				"t.sequence=s.name AND t.id IN (%s)"
			)
			sql = sql % (self.table, ",".join(map(str, self.ids.values())))

		current = 0
		seqs = None

		insert_sql = "INSERT INTO primer VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
		
		self.db.begin()
		for item in self.db.query(sql):
			if seqs is None or item[2] not in seqs:
				seqs = fasta.GzipFasta(item[0])
			
			start = item[3] - self.flank
			if start < 1:
				start = 1
			end = item[4] + self.flank
			
			target = dict(
				SEQUENCE_ID = "%s-%s" % (self.table, item[1]),
				SEQUENCE_TEMPLATE = seqs.get_seq_by_loci(item[2], start, end),
				SEQUENCE_TARGET = [item[3]-start, item[5]],
				SEQUENCE_INTERNAL_EXCLUDED_REGION = [item[3]-start, item[5]]
			)

			primerdesign.setSeqArgs(target)
			res = primerdesign.runDesign(False)
			current += 1

			primer_count = res['PRIMER_PAIR_NUM_RETURNED']
			for i in range(primer_count):
				primer = [None, target['SEQUENCE_ID'], i+1]
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
		
			self.emit_progress(int(current/total_ids*100))

		self.db.commit()

		self.emit_finish('Primer design completed')

class ExportFastaWorker(Worker):
	def __init__(self, table, ids, flank, outfile):
		super(ExportFastaWorker, self).__init__()
		self.table = table
		self.ids = ids
		self.flank = flank
		self.outfile = outfile

	def process(self):
		total_ssrs = self.db.get_one("SELECT COUNT(1) FROM %s" % self.table)
		total_ids = len(self.ids)

		if total_ssrs == total_ids:
			sql = (
				"SELECT t.*, f.path FROM %s AS t,fasta AS f,seq AS s "
				"WHERE f.id=s.fid AND t.sequence=s.name"
			)
			sql = sql % self.table
		else:
			sql = (
				"SELECT t.*, f.path FROM %s AS t,fasta AS f,seq AS s "
				"WHERE f.id=s.fid AND t.sequence=s.name AND t.id IN (%s)"
			)
			sql = sql % (self.table, ",".join(map(str, self.ids.values())))

		current = 0
		seqs = None

		fp = open(self.outfile, 'w')

		for item in self.db.query(sql):
			if seqs is None or item.sequence not in seqs:
				seqs = fasta.GzipFasta(item.path)
			
			start = item.start - self.flank
			if start < 1:
				start = 1
			end = item.end + self.flank
			ssr = seqs.get_seq_by_loci(item.sequence, start, end)
			name = ">%s%s %s:%s-%s|motif:%s" % (self.table.upper(), item.id, item.sequence, item.start, item.end, item.motif)
			fp.write("%s\n%s" % (name, format_fasta_sequence(ssr, 70)))
			
			current += 1
			self.emit_progress(int(current/total_ids*100))

		fp.close()

		self.emit_finish("Export fasta completed.")

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
		self.emit_message("Building interval tree")
		interval_forest = {}
		genes_info = {}
		f = check_gene_annot_format(self.annot_file)
		if f == 'GFF':
			features = get_gff_coordinate(self.annot_file)
		else:
			features = get_gtf_coordinate(self.annot_file)

		self.emit_message("Building interval tree")
		for feature in features:
			if feature[0] == 'GENE':
				genes_info[feature[1]] = feature[2]
				continue

			if feature[1] not in interval_forest:
				interval_forest[feature[1]] = intersection.IntervalTree()

			interval_forest[feature[1]].insert(feature[2], feature[3], (feature[0], feature[4]))

		total = self.db.get_one("SELECT COUNT(1) FROM %s LIMIT 1" % self.table)
		current = 0
		for ssr in self.db.get_cursor().execute("SELECT * FROM %s" % self.table):
			self.emit_message("Locating %ss in sequence %s" % (self.table.upper(), ssr.sequence))
			current += 1
			progress = int(current/total*100)
			if progress%2 == 0:
				self.emit_progress(progress)

			if ssr.sequence not in interval_forest:
				continue

			res = interval_forest[ssr.sequence].find(ssr.start, ssr.end)
			if not res: continue

			res = {f:g for f, g in res}

			record = [None, self.table, ssr.id]
			for feat in ['CDS', '5UTR', '3UTR', 'EXON', 'INTRON']:
				if feat in res:
					gid = res[feat]
					record.extend([gid, genes_info[gid], feat])
					self.db.get_cursor().execute("INSERT INTO location VALUES (?,?,?,?,?,?)", record)
					break

		self.emit_message("Creating query index")
		self.db.query("CREATE INDEX loci ON location (target, category)")
		self.db.query("CREATE INDEX sel ON location (category, feature)")

		self.emit_finish("%s location completed." % self.table)

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
		message = "Download %s fasta completed" % self.acc
		try:
			r = requests.get(url, timeout=10, stream=True)
		except:
			message = "Can not connect to NCBI server"
		else:
			if r.status_code == requests.codes.ok:
				self.total = 0
				self.start = time.time()
				with open(self.outfile, "wb") as fh:
					try:
						for chunk in r.iter_content(chunk_size=1024):
							self.total += len(chunk)
							fh.write(chunk)
							self.emit_message(self.progressing())
					except Exception, e:
						message = str(e)
			else:
				message = "%s error: %s" % (r.status_code, r.reason)

		self.emit_finish(message)

	def progressing(self):
		total = human_size(self.total)
		timer = time.time() - self.start
		speed = human_size(self.total/timer)
		return "Downloaded %s, Speed %s/s" % (total, speed)
