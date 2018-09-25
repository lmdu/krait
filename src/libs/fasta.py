#!/usr/bin/env python
import os
import gzip
import mmap
import collections
from . import kseq
from . import indexed_gzip as igzip

class GzipFasta:
	def __init__(self, fasta_file, rebuild=False):
		self.fasta_file = fasta_file
		self.index_file = "%s.fidx" % fasta_file
		self.rebuild = rebuild
		self.buff = {'name': None, 'seq': None}
		
		self._read_index()
		
		#open fasta file
		kseq.open_fasta(self.fasta_file)
	
	def __iter__(self):		
		return self

	def __next__(self):
		seq = kseq.iter_seq()
		if seq is None:
			#end loop and close fasta
			kseq.close_fasta()
			raise StopIteration
		return seq

	def __len__(self):
		return len(self._index)

	def __contains__(self, name):
		return name in self._index

	def __getitem__(self, name):
		return self.get_seq_by_name(name)

	def _read_fasta(self):
		if self.fasta_file.endswith('.gz'):
			idxfile = '{}.gzidx'.format(self.fasta_file)
			if os.path.exists(idxfile):
				handler = igzip.IndexedGzipFile(self.fasta_file, index_file=idxfile)
			else:
				handler = igzip.IndexedGzipFile(self.fasta_file)
				#handler.build_full_index()
				#handler.export_index(idxfile)
		else:
			f = open(self.fasta_file, 'rb')
			handler = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
		
		return handler
	
	def _read_index(self):
		self._index = collections.OrderedDict()
		if self.rebuild or not os.path.exists(self.index_file) or not os.path.getsize(self.index_file):
			self._build_index()

		with open(self.index_file) as fh:
			for line in fh:
				cols = line.strip().split('\t')
				self._index[cols[0]] = [int(x) for x in cols[1:]]

	def _build_index(self):
		chroms = kseq.build_index(self.fasta_file)
		with open(self.index_file, 'w') as fw:
			for chrom in chroms:
				fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % chrom)

	def keys(self):
		for k in self._index:
			yield k

	def get_gc(self, name):
		return self._index[name][3]

	def get_ns(self, name):
		return self._index[name][4]

	def get_len(self, name):
		'''
		get the sequence length by the name
		@para name str, the sequence name
		@return int, length
		'''
		return self._index[name][2]

	def get_total_length(self):
		return sum(self.get_len(name) for name in self._index)

	def get_seq_by_name(self, name):
		'''
		get whole seqence from fasta file by the name
		@para name str, the fasta sequence name
		@return str, the dna sequence
		'''
		start, length = self._index[name][:2]
		fp = self._read_fasta()
		fp.seek(start)
		content = fp.read(length)
		fp.close()
		return kseq.clean_seq(content.decode('utf-8'))

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

	def get_seq_by_pos(self, name, start, end):
		start, length = self._index[name][:2]
		fp = self._read_fasta()
		fp.seek(start)
		content = fp.read(length)
		fp.close()
		return kseq.sub_seq(content.decode('utf-8'), start, end)


if __name__ == '__main__':
	GzipFasta('D:/tandem/data/Danio_rerio.GRCz10.dna.chromosome.1.fa', True)
