#!/usr/bin/env python
import gzip
import apsw
import operator

class Data(dict):
	def __getattr__(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

	def __setattr__(self, name, val):
		self[name] = val

create_table_sql = '''
CREATE TABLE parent (
	id integer primary key,
	seqid text,
	type text,
	start integer,
	end integer,
	strand text,
	featid text,
	name text
);
CREATE TABLE parent_relation (
	child integer,
	parent integer
);
CREATE TABLE feature (
	id integer primary key,
	seqid text,
	type text,
	start integer,
	end integer,
	strand text,
	featid text,
	name text
);
CREATE TABLE feat_relation (
	feature integer,
	parent integer
);
'''

db = apsw.Connection(':memory:')
db.cursor().execute("pragma synchronous=0;")
db.cursor().execute(create_table_sql)

class GFFParser(object):
	def __init__(self, handler, _format):
		self.handler = handler
		self.format = _format

	def _parse_attrs(self, attrstr):
		attrs = Data()
		if self.format == 'GFF':
			for item in attrstr.split(';'):
				name, value = item.split('=')
				attrs[name] = value
		else:
			for item in attrstr.split(';'):
				if not item:
					continue
				name, value = item.strip().strip('"').split('"')
				attrs[name.strip()] = value

		return attrs

	def __iter__(self):
		while 1:
			self.line = self.handler.readline()
			if not self.line:
				break
			
			if self.line[0] == '#':
				continue
			else:
				break

		return self

	def __next__(self):
		if not self.line:
			raise StopIteration
		
		cols = self.line.strip().split('\t')
		row = Data()
		row.seqid = cols[0]
		row.feature = cols[2].lower()
		row.start = int(cols[3])
		row.end = int(cols[4])
		row.strand = cols[6]
		row.attrs = self._parse_attrs(cols[-1])
		
		#read a new line
		self.line = self.handler.readline()

		return row


class Annotation(object):
	def __init__(self, annot_file):
		self.annot_file = annot_file
		self.format = self._check_format()

	def _open(self, infile):
		if self.annot_file.endswith('.gz'):
			return gzip.open(infile, 'rt')
		else:
			return open(infile)

	def _check_format(self):
		with self._open(self.annot_file) as fh:
			for line in fh:
				if line[0] == '#':
					continue

				cols = line.strip().split('\t')

				if len(cols) == 9:
					if cols[-1].count('=') >= 2:
						return 'GFF'

					elif 'gene_id' in cols[-1]:
						return 'GTF'
					
					else:
						raise Exception('The annotation file is not GTF or GFF format')

				else:
					raise Exception('The annotation file is not GTF or GFF format')

	#@profile
	'''
	def parse_gff(self):
		prev_parent = None
		exons = []
		with self._open(self.annot_file) as fh:
			for row in GFFParser(fh, self.format):
				if row.feature == 'region':
					continue

				if 'Name' not in row.attrs:
					row.attrs.Name = row.attrs.ID

				feature = Feature.create(
					seqid = row.seqid,
					type = row.feature,
					start = row.start,
					end = row.end,
					strand = row.strand,
					featid = row.attrs.ID,
					name = row.attrs.Name
				)

				if 'Parent' in row.attrs:
					try:
						parent = Feature.get(Feature.featid==row.attrs.Parent)
					except Feature.DoesNotExist:
						raise Exception('Can not find parent for {} (ID:{}) in annotation file'.format(row.feature, row.attrs.ID))

					Relation.create(feature=feature, parent=parent)

				if row.feature == 'exon':
					if prev_parent == row.attrs.Parent:
						exons.append(row)
					else:
						if exons:
							self.get_introns(exons, parent)
						
						exons = [row]
						prev_parent = row.attrs.Parent

			if exons:
				self.get_introns(exons, parent)
	'''
	@profile
	def parse_gtf(self):
		prev_parent = None
		exons = []
		genes = {}

		with self._open(self.annot_file) as fh:
			for row in GFFParser(fh, self.format):
				if row.feature == 'gene':
					if 'gene_name' not in row.attrs:
						row.attrs['gene_name'] = row.attrs.gene_id
					
					feature = (None, row.seqid, row.feature, row.start, row.end, row.strand, row.attrs.gene_id, row.attrs.gene_name)

					db.cursor().execute("INSERT INTO parent VALUES (?,?,?,?,?,?,?,?)", feature)

					genes[row.attrs.gene_id] = db.last_insert_rowid()

				elif row.feature == 'transcript':
					pass

				else:
					feature = (None, row.seqid, row.feature, row.start, row.end, row.strand, '', '')
					db.cursor().execute("INSERT INTO feature VALUES (?,?,?,?,?,?,?,?)", feature)
					db.cursor().execute("INSERT INTO feat_relation VALUES (?,?)", (db.last_insert_rowid(), genes[row.attrs.gene_id]))
					
					if row.feature == 'exon':
						if prev_parent == row.attrs.transcript_id:
							exons.append(row)
						else:
							if exons:
								self.get_introns(exons, genes)
							
							exons = [row]
							prev_parent = row.attrs.transcript_id

			if exons:
				self.get_introns(exons, genes)


	def get_introns(self, exons, genes):
		exons = sorted(exons, key=operator.attrgetter('start'))
		parent = genes[exons[0].attrs.gene_id]
		for idx, exon in enumerate(exons):
			start, end = exon.start, exon.end

			if idx < len(exons)-1:
				start = end + 1
				end = exons[idx+1].start - 1
				feature = (None, exons[0].seqid, 'intron', start, end, exons[0].strand, '', '')
				db.cursor().execute("INSERT INTO feature VALUES (?,?,?,?,?,?,?,?)", feature)
				db.cursor().execute("INSERT INTO feat_relation VALUES (?,?)", (db.last_insert_rowid(), genes[exons[0].attrs.gene_id]))


	def parse_annot(self):
		if self.format == 'GTF':
			self.parse_gtf()
		else:
			self.parse_gff()

	'''
	def get_parent(self, feature):
		relate1 = Relation.get(Relation.feature==feature)
		try:
			relate2 = Relation.get(Relation.feature==relate1.parent)
		except Relation.DoesNotExist:
			return (relate1.parent.featid, relate1,parent.name)
		else:
			return (relate2.parent.featid, relate2.parent.name)

	def get_feature(self, chrom, start, end):
		features = set()
		for feature in Feature.select().where(Feature.seqid==chrom, Feature.start<=start, Feature.end>=end):
			if feature.type in ('gene', 'mrna', 'transcript'):
				continue
			fid, name = self.get_parent(feature)

			features.add(feature.type, fid, name)
	'''

if __name__ == '__main__':
	annot = Annotation('/mnt/d/research/tandem/data/Danio_rerio.GRCz10.89.gtf.gz')
	annot.parse_annot()

	#annot.get_feature('4', 48623, 48670)

