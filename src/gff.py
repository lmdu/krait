import gzip
import numpy

from libs import ncls
from utils import Data

def check_gene_annot_format(annot_file):
	if annot_file.endswith('.gz'):
		fh = gzip.open(annot_file, 'rt')
	else:
		fh = open(annot_file)

	for line in fh:
		if line[0] == '#':
			continue

		cols = line.strip().split('\t')

		if len(cols) != 9:
			raise Exception('The annotation file is not GTF or GFF format')
		
		if cols[-1].count('=') >= 1:
			return 'GFF'

		elif 'gene_id "' in cols[-1]:
			return 'GTF'

		else:
			raise Exception('The annotation file is not GTF or GFF format')

class AnnotParser(object):
	def __init__(self, annot_file):
		self.annot_file = annot_file
		self.gene_mapping = {}
		self.gene_info = []
		self.feature_mapping = {}
		self.interval_forest = {}

		self.repeat_types = {'ssr': 1, 'cssr': 2, 'issr': 3, 'vntr': 4}
		self.featid_mapping = {'CDS': 1, 'exon': 2, '3UTR': 3, 'intron': 4, '5UTR': 5}
		
		self.get_gene_mapping()
		self.create_interval_tree()

	def __contains__(self, item):
		return item in self.interval_forest

	def parse(self):
		if self.annot_file.endswith('.gz'):
			fh = gzip.open(self.annot_file, 'rt')
		else:
			fh = open(self.annot_file)

		for line in fh:
			if line[0] == '#': continue
			
			cols = line.strip().split('\t')
			
			record = Data(
				seqid = cols[0],
				feature = cols[2].upper(),
				start = int(cols[3]),
				end = int(cols[4]),
				attrs = Data()
			)
			
			for item in cols[-1].split(';'):
				if not item:
					continue
				
				#if _format == 'GFF':
				#	name, value = item.split('=')
				#else:
				#	name, value = item.strip().strip('"').split('"')

				name, value = self.split_val(item)
				
				record.attrs[name.strip().upper()] = value
			
			yield record

		fh.close()

	def split_value(self, item):
		pass

	def get_gene_mapping(self):
		pass

	def get_features(self):
		pass

	def create_interval_tree(self):
		prev_chrom = None
		starts = []
		ends = []
		indexes = []
		feat_id = 0

		for feature in self.get_features():
			feat_id += 1
			self.feature_mapping[feat_id] = feature[3:]

			if feature[0] != prev_chrom:
				if starts:
					starts = numpy.array(starts, dtype=numpy.long)
					ends = numpy.array(ends, dtype=numpy.long)
					indexes = numpy.array(indexes, dtype=numpy.long)
					self.interval_forest[prev_chrom] = ncls.NCLS(starts, ends, indexes)

				prev_chrom = feature[0]
				starts = []
				ends = []
				indexes = []

			starts.append(feature[1])
			ends.append(feature[2])
			indexes.append(feat_id)

		if starts:
			starts = numpy.array(starts, dtype=numpy.long)
			ends = numpy.array(ends, dtype=numpy.long)
			indexes = numpy.array(indexes, dtype=numpy.long)
			self.interval_forest[prev_chrom] = ncls.NCLS(starts, ends, indexes)

	def mapping(self, chrom, start, end):
		if chrom not in self.interval_forest:
			return None

		res = set(self.interval_forest[chrom].find_overlap(start, end))

		if not res:
			return None

		feats = [self.feature_mapping[feat_id[2]] for feat_id in res]

		for candidate in ['CDS', 'exon', 'UTR', 'intron']:
			for feat, gid in feats:
				if candidate in feat:
					return (self.featid_mapping[feat], self.gene_mapping[gid])

		return None

class GFFParser(AnnotParser):
	def split_val(self, item):
		return item.split('=')

	def get_gene_mapping(self):
		gene_num = 0
		for row in self.parse():
			if row.feature == 'REGION':
				continue

			if row.feature != 'GENE':
				if 'PARENT' in row.attrs:
					continue

			if 'ID' in row.attrs:
				gene_id = row.attrs.ID

			elif 'GENE' in row.attrs:
				gene_id = row.attrs.GENE

			elif 'NAME' in row.attrs:
				gene_id = row.attrs.NAME

			else:
				raise Exception(row)

			if 'NAME' in row.attrs:
				gene_name = row.attrs.NAME

			elif 'PRODUCT' in row.attrs:
				gene_name = row.attrs.PRODUCT

			elif 'GENE' in row.attrs:
				gene_name = row.attrs.GENE

			elif 'ID' in row.attrs:
				gene_name = row.attrs.ID

			else:
				raise Exception(row)

			biotype = row.attrs.get('GENE_BIOTYPE', row.feature)

			gene_num += 1

			self.gene_mapping[gene_id] = gene_num

			self.gene_info.append((gene_num, row.seqid, row.start, row.end, gene_id, gene_name, biotype))

	def get_features(self):
		father = None
		exons = []

		parents = {}

		for r in self.parse():
			if r.feature == 'REGION':
				continue

			elif r.feature == 'GENE':
				if 'ID' in r.attrs:
					parents[r.attrs.ID] = r.attrs.ID
				elif 'GENE' in r.attrs:
					parents[r.attrs.GENE] = r.attrs.GENE
					parents['gene-{}'.format(r.attrs.GENE)] = r.attrs.GENE
				elif 'NAME' in r.attrs:
					parents[r.attrs.NAME] = r.attrs.NAME

			elif r.feature == 'CDS':
				if 'PARENT' in r.attrs:
					yield (r.seqid, r.start, r.end, 'CDS', parents[r.attrs.PARENT])
				else:
					yield (r.seqid, r.start, r.end, 'CDS', r.attrs.ID)
			
			elif r.feature == 'FIVE_PRIME_UTR':
				if 'PARENT' in r.attrs:
					yield (r.seqid, r.start, r.end, '5UTR', parents[r.attrs.PARENT])
				else:
					yield (r.seqid, r.start, r.end, '5UTR', r.attrs.ID)

			elif r.feature == 'THREE_PRIME_UTR':
				if 'PARENT' in r.attrs:
					yield (r.seqid, r.start, r.end, '3UTR', parents[r.attrs.PARENT])
				else:
					yield (r.seqid, r.start, r.end, '3UTR', r.attrs.ID)
			
			elif r.feature == 'EXON':
				try:
					mother = r.attrs.PARENT
				except AttributeError:
					continue

				if father == mother:
					exons.append((r.seqid, r.start, r.end, 'exon', parents[r.attrs.PARENT]))
				
				else:
					if exons:
						exons = sorted(exons, key=lambda x: x[2])
						for idx, exon in enumerate(exons):
							yield exon

							if idx < len(exons)-1:
								start = exon[2] + 1
								end = exons[idx+1][1] - 1
								yield (exons[0][0], start, end, 'intron', parents[r.attrs.PARENT])
					
					exons = [(r.seqid, r.start, r.end, 'exon', parents[r.attrs.PARENT])]
					father = mother
			
			else:
				if 'ID' in r.attrs:
					try:
						parents[r.attrs.ID] = parents[r.attrs.PARENT]
					except:
						parents[r.attrs.ID] = r.attrs.ID

		exons = sorted(exons, key=lambda x: x[2])
		
		for idx, exon in enumerate(exons):
			yield exon

			if idx < len(exons)-1:
				start = exon[2] + 1
				end = exons[idx+1][1] - 1
				yield (exons[0][0], start, end, 'intron', exons[0][4])

class GTFParser(AnnotParser):
	def split_val(self, item):
		return item.strip().strip('"').split('"')

	def get_gene_mapping(self):
		gene_num = 0
		for row in self.parse():
			if row.feature is not 'GENE':
				continue

			gene_id = row.attrs.GENE_ID
			gene_name = row.attrs.get('GENE_NAME', row.attrs.GENE_ID)
			biotype = row.attrs.get('GENE_BIOTYPE', row.feature)

			gene_num += 1
			self.gene_mapping[gene_id] = gene_num

			self.gene_info.append((gene_num, row.seqid, row.start, row.end, gene_id, gene_name, biotype))

	def get_features(self):
		father = None
		exons = []
		for row in self.parse():
			parent = row.attrs.GENE_ID

			if row.feature == 'CDS':
				yield (r.seqid, r.start, r.end, 'CDS', parent)
			
			elif r.feature == 'FIVE_PRIME_UTR':
				yield (r.seqid, r.start, r.end, '5UTR', parent)
			
			elif r.feature == 'THREE_PRIME_UTR':
				yield (r.seqid, r.start, r.end, '3UTR', parent)
			
			elif r.feature == 'UTR':
				yield (r.seqid, r.start, r.end, 'UTR', parent)
			
			elif r.feature == 'EXON':
				mother = r.attrs.TRANSCRIPT_ID

				if father == mother:
					exons.append((r.seqid, r.start, r.end, 'exon', parent))

				else:
					if exons:
						exons = sorted(exons, key=lambda x: x[1])
						
						for idx, exon in enumerate(exons):
							yield exon

							if idx < len(exons)-1:
								start = exon[2] + 1
								end = exons[idx+1][1] - 1
								yield (exons[0][0], start, end, 'intron', exons[0][4])
					
					exons = [(r.seqid, r.start, r.end, 'exon', parent)]
					father = mother

		if exons:
			exons = sorted(exons, key=lambda x: x[1])

			for idx, exon in enumerate(exons):
				yield exon

				if idx < len(exons)-1:
					start = exon[2] + 1
					end = exons[idx+1][1] - 1
					yield (exons[0][0], start, end, 'intron', exons[0][4])


