#!/usr/bin/env python
import os
import gzip
import collections

class Fasta:
	def __init__(self, fasta_file, rebuild=False):
		self.fasta_file = fasta_file
		self.index_file = "%s.idx" % fasta_file
		self.rebuild = rebuild
		self._check_fasta()
		self._read_index()

	def __iter__(self):
		for name in self._index:
			seq = self.get_sequence(name)
			yield (name, seq)

	def __len__(self):
		return len(self._index)

	def _read_fasta(self):
		if self.fasta_file.endswith('.gz'):
			return gzip.open(self.fasta_file, 'rb')
		else:
			return open(self.fasta_file, 'rb')

	def _check_fasta(self):
		fh = self._read_fasta()
		line = fh.readline()
		if line[0] != '>':
			raise Exception("%s is not fasta format file" % self.fasta_file)
		fh.close()
	
	def _read_index(self):
		self._index = collections.OrderedDict()
		if not os.path.exists(self.index_file) or self.rebuild:
			self._build_index()

		with open(self.index_file) as fh:
			for line in fh:
				cols = line.strip().split('\t')
				self._index[cols[0]] = (int(cols[1]), int(cols[2]))

	def _build_index(self):
		fh = self._read_fasta()
		line = fh.readline()
		name = line[1:].strip().split()[0]
		start = fh.tell()

		with open(self.index_file, 'w') as fw:
			while 1:
				line = fh.readline()
				if not line:
					length = fh.tell() - start
					fw.write("%s\t%s\t%s\n" % (name, start, length))
					break

				if line[0] == '>':
					pos = fh.tell()
					length = pos - start - len(line)
					fw.write("%s\t%s\t%s\n" % (name, start, length))
					start =  fh.tell()
					name = line[1:].strip().split()[0]
		fh.close()

	def keys(self):
		return self._index.keys()

	def get_sequence(self, name):
		start, length = self._index[name]
		fp = self._read_fasta()
		fp.seek(start)
		content = fp.read(length)
		fp.close()
		return content.replace('\n', '').replace('\r', '').upper()


