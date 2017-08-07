#!/usr/bin/env python
import os
import gzip
import kseq
import collections

class GzipFasta:
	def __init__(self, fasta_file, rebuild=False):
		self.fasta_file = fasta_file
		self.index_file = "%s.idx" % fasta_file
		self.rebuild = rebuild
		self.buff = {'name': None, 'seq': None}
		self._read_index()

	def __enter__(self):
		kseq.open_fasta(self.fasta_file)
		return self

	def __exit__(self, exc_type, exc_value, exc_tb):
		if exc_tb is None:
			kseq.close_fasta()
		else:
			raise Exception("Close fasta error")

	def __iter__(self):
		return self

	def next(self):
		seq = kseq.iter_seq()
		if seq is None:
			raise StopIteration
		else:
			return seq

	def __len__(self):
		return len(self._index)

	def __contains__(self, name):
		return name in self._index

	def __getitem__(self, name):
		return self.get_seq_by_name(name)

	def _read_fasta(self):
		if self.fasta_file.endswith('.gz'):
			return gzip.open(self.fasta_file, 'rb')
		else:
			return open(self.fasta_file, 'rb')
	
	def _read_index(self):
		self._index = collections.OrderedDict()
		if not os.path.exists(self.index_file) or self.rebuild:
			self._build_index()

		with open(self.index_file) as fh:
			for line in fh:
				cols = line.strip().split('\t')
				self._index[cols[0]] = (int(cols[1]), int(cols[2]))

	def _build_index(self):
		chroms = kseq.build_index(self.fasta_file)
		with open(self.index_file, 'w') as fw:
			for chrom in chroms:
				fw.write("%s\t%s\t%s\n" % chrom)

	@property
	def keys(self):
		return self._index.keys()

	def get_seq_by_name(self, name):
		'''
		get whole seqence from fasta file by the name
		@para name str, the fasta sequence name
		@return str, the dna sequence
		'''
		start, length = self._index[name]
		fp = self._read_fasta()
		fp.seek(start)
		content = fp.read(length)
		fp.close()
		return kseq.clean_seq(content)

	def get_seq_by_loci(self, name, start, end):
		'''
		get sequence fragment from a fasta sequence by start and end position
		@para name str, the fasta name
		@para start int, the start position
		@para end int, the end position
		@return str, the extracted sequence
		'''
		if self.buff['name'] != name:
			self.buff['name'] = name
			self.buff['seq'] = self.get_seq_by_name(name)

		return self.buff['seq'][start-1:end]

if __name__ == '__main__':
	GzipFasta('D:/tandem/data/Danio_rerio.GRCz10.dna.chromosome.1.fa', True)
