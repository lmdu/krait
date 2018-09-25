import os
import sys
import csv
import glob
import time
import motif
import random
import argparse
import multiprocessing as mp

from libs import *

class Workers(object):
	def __init__(self, cpus):
		#CPU counts specified by user
		self.cpus = cpus
		
		#store results
		self.results = {}
		
		#number of current tasks in pool
		self.tasks = 0

		if self.cpus > mp.cpu_count():
			self.cpus = mp.cpu_count()

		self.pool = mp.Pool(self.cpus)

	def add_task(self, job, args):
		self.pool.apply_async(
			func = job,
			args = args,
			callback = self.success,
			error_callback = self.failure
		)
		self.tasks += 1

	def success(self, res):
		self.tasks -= 1
		if res:
			self.results[res[0][0]] = res

	def failure(self, error):
		print("Error occured: {}".format(error))

	def full(self):
		return self.tasks >= self.cpus

	def release(self):
		self.pool.close()
		self.pool.join()

def format_to_gff(feature, row):
	cols = [row.sequence, 'Krait', feature, row.start, row.end, '.', '+', '.', []]
	cols[-1].append("ID={}{}".format(feature, row.id))
	cols[-1].append("Motif={}".format(row.motif))
	for k in row.getKeys():
		if k not in ['id', 'sequence', 'start', 'end', 'motif']:
			cols[-1].append("{}={}".format(k.capitalize(), row.value(k)))
	cols[-1] = ";".join(cols[-1])
	return cols

def save_result(fw, outfmt, rows, ssr_type='SSR'):
	if outfmt == 'csv':
		writer = csv.writer(fw)
	else:
		writer = csv.writer(fw, delimiter='\t')

	if outfmt == 'gff':
		write_line = lambda x: writer.writerow(format_to_gff(ssr_type, x))
	else:
		write_line = lambda x: writer.writerow(x)

	for row in rows:
		write_line(row)


def search_ssr(name, seq, repeats, level):
	motifs = motif.StandardMotif(level)
	ssrs = tandem.search_ssr(seq, repeats)
	return [[name, motifs.standard(ssr[0])] + list(ssr) for ssr in ssrs]

def concatenate_cssr(cssrs):
	seqname = cssrs[-1][0]
	start = cssrs[0][5]
	end = cssrs[-1][6]
	complexity = len(cssrs)
	motif = "-".join([cssr[2] for cssr in cssrs])
	length = sum(cssr[7] for cssr in cssrs)
	gap = sum(cssr[5]-cssrs[idx][6]-1 for idx, cssr in enumerate(cssrs[1:]))
	structure = "-".join(["(%s)%s" % (cssr[2], cssr[4]) for cssr in cssrs])
	return (seqname, start, end, motif, complexity, length, gap, structure)

def search_cssr(name, seq, repeats, dmax):
	ssrs = tandem.search_ssr(seq, repeats)
	res = []
	cssrs = [ssrs[0]]
	for ssr in ssrs[1:]:
		d = ssr[5] - cssrs[-1][6] - 1
		if d <= dmax:
			cssrs.append(ssr)
		else:
			if len(cssrs) > 1:
				res.append(concatenate_cssr(cssrs))

			cssrs = [ssr]

	if len(cssrs) > 1:
		res.append(concatenate_cssr(cssrs))

	return res

def search_issr(name, seq):
	pass

def search_vntr(name, seq):
	pass


class Jobs(object):
	def __init__(self, args):
		self.args = args

		#open fasta file
		self.seqs = fasta.GzipFasta(args.infile)

		#create multiple process pool
		self.pool = Workers(self.args.cpus)

		#start process job
		self.run_jobs()

		#result row number
		self.row_num = 0

		#save result to file or terminal
		self.save_result()

	def __getstate__(self):
		self_dict = self.__dict__.copy()
		del self_dict['pool']
		return self_dict

	def __setstate__(self, state):
		self.__dict__.update(state)

	def run_jobs(self):
		target = self.get_func()
		while 1:
			if not self.pool.full():
				job = self.get_job()

				if job is None:
					break
				
				self.pool.add_task(target, job)
			else:
				time.sleep(0.01)

		self.pool.release()

		self.save_result()
	
	def save_result(self):
		if self.args.outfile == 'stdout':
			fw = sys.stdout
		else:
			fw = open(self.args.outfile, 'w', newline='')

		if self.args.outfmt == 'csv':
			writer = csv.writer(fw)
		else:
			writer = csv.writer(fw, delimiter='\t')

		if self.args.outfmt == 'gff':
			write_line = lambda x: writer.writerow(format_gff(x))
		else:
			write_line = lambda x: writer.writerow(x)

		for k in self.seqs.keys():
			if k in self.pool.results:
				for row in self.pool.results[k]:
					write_line(row)

		if self.args.outfile == 'stdout':
			fw.flush()
		else:
			fw.close()

	def get_func(self):
		pass

	def get_job(self):
		pass

	def format_gff(self, row):
		pass

class SSRSearchJob(Jobs):
	def __init__(self, args):
		super(SSRSearchJob, self).__init__(args)

	def get_func(self):
		return search_ssr

	def get_job(self):
		try:
			name, seq = next(self.seqs)
		except StopIteration:
			return None
		
		return (name, seq, self.args.repeats, self.args.level)

	def format_gff(self, row):
		types = {1:'Mono', 2:'Di', 3:'Tri', 4:'Tetra', 5:'Penta', 6:'Hexa'}
		self.row_num += 1
		attrs = 'ID={};Motif={};Standard={};Type={};Repeat={};Length={}'.format(
			self.row_num, row[2], row[1], types[row[3]], row[4], row[7])
		fields = [row[0], 'Krait', self.args.ssr_type.upper(), row[5], row[6], '.', '.', '.', attrs]
		return fields

class CSSRSearchJob(Jobs):
	def __init__(self, args):
		super(CSSRSearchJob, self).__init__(args)

	def get_func(self):
		return search_cssr

	def get_job(self):
		try:
			name, seq = next(self.seqs)
		except StopIteration:
			return None
		
		return (name, seq, self.args.repeats, self.args.dmax)

	def format_gff(self, row):
		self.row_num += 1
		attrs = 'ID={};Motif={};Complexity={};Length={};Gap={};Structure={}'.format(
			self.row_num, row[3], row[4], row[5], row[6], row[7])
		fields = [row[0], 'Krait', self.args.ssr_type.upper(), row[1], row[2], '.', '.', '.', attrs]
		return fields

def search_tandem(args):
	if args.ssr_type == 'ssr':
		func = search_ssr
	
	elif args.ssr_type == 'cssr':
		func = search_cssr

	seqs = fasta.GzipFasta(args.infile)
	items = search_parameters(seqs, args)

	worker = Workers(func, items, args.cpus)
	
	if outfile == 'stdout':
		fw = sys.stdout
	else:
		fw = open(outfile, 'w', newline='')

	for k in seqs.keys():
		if k in worker.results:
			save_result(fw, args.outfmt, worker.results[k], args.ssr_type.upper())

	if args.outfile == 'stdout':
		fw.flush()
	else:
		fw.close()


def extract_flank_sequence():
	pass

def mapping_to_gene():
	pass

def design_primer():
	pass



if __name__ == '__main__':
	mp.freeze_support()

	description = '''Search for microsatellite from large genome

Cite:
    Du L, Zhang C, Liu Q, Zhang X, Yue B (2018) Krait: an
    ultrafast tool for genome-wide survey of microsatellites
    and primer design. Bioinformatics, 34(4):618-683.
	'''

	parser = argparse.ArgumentParser(
		prog = 'krait',
		usage = 'krait COMMAND [OPTIONS]',
		description = description,
		epilog = '''Contact:
	    Lianming Du (dulianming@cdu.edu.cn)
		''',
		formatter_class = argparse.RawDescriptionHelpFormatter
	)

	parser.add_argument('-v', '--version',
		action = 'version',
		version = '%(prog)s version 1.0'
	)

	subparsers = parser.add_subparsers(
		title = 'Commands',
		prog = 'krait',
		metavar = '',
	)
	parser_search = subparsers.add_parser('search',
		help = "Search for simple tandem repeats",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_search.set_defaults(cmd='search')

	parser_search.add_argument('-t', '--type',
		dest = 'ssr_type',
		choices = ('ssr', 'cssr', 'issr', 'vntr'),
		default = 'ssr',
		metavar = '',
		help = 'ssr type (ssr, cssr, issr, vntr)'
	)
	parser_search.add_argument('-i', '--in',
		dest = 'infile',
		metavar = 'fasta',
		required = False,
		type = lambda x: os.path.normcase(x),
		help = "input fasta file or gzipped fasta file"
	)
	parser_search.add_argument('-o', '--out',
		dest = 'outfile',
		default = 'stdout',
		metavar = '',
		type = lambda x: os.path.normcase(x),
		help = "output file name"
	)
	parser_search.add_argument('-f', '--format',
		dest = 'outfmt',
		default = 'tsv',
		choices = ('tsv', 'csv', 'gff'),
		metavar = '',
		help="output format [tsv, csv, gff]"
	)
	parser_search.add_argument('-c', '--cpus',
		dest = 'cpus',
		default = 1,
		type = int,
		metavar = '',
		help = 'number of threads'
	)
	parser_search.add_argument('-l', '--level',
		dest = 'level',
		default = 3,
		type = int,
		metavar = '',
		help = 'motif standardization level'
	)

	ssr_group = parser_search.add_argument_group('SSR', 'perfect SSR search parameter')
	ssr_group.add_argument('--repeats',
		dest = 'repeats',
		nargs = 6,
		default = [12,7,5,4,4,4],
		metavar = ('mono', 'di', 'tri', 'tetra', 'penta', 'hexa'),
		type = int,
		help = "minimum repeats for mono di tri tetra penta hexa SSR search"
	)

	cssr_group = parser_search.add_argument_group('cSSR', 'compound SSR search parameter')
	cssr_group.add_argument('--dmax',
		dest = 'dmax',
		default = 10,
		metavar = '',
		type = int,
		help = "maximum distance between adjacent perfect SSRs"
	)

	issr_group = parser_search.add_argument_group('iSSR', 'imperfect SSR search parameter')
	issr_group.add_argument('--min-seed-rep',
		dest = 'min_seed_repeats',
		default = 3,
		metavar = '',
		type = int,
		help = "minimum seed repeats"
	)
	issr_group.add_argument('--min-seed-len',
		dest = 'min_seed_length',
		default = 8,
		metavar = '',
		type = int,
		help = "minimum seed length"
	)
	issr_group.add_argument('--max-edits',
		dest = 'max_consecutive_edits',
		default = 3,
		metavar = '',
		type = int,
		help = "maximum consecutive edits"
	)
	issr_group.add_argument('--mis-penalty',
		dest = 'mis_penalty',
		default = 1,
		metavar = '',
		type = int,
		help = "mismatch penalty"
	)
	issr_group.add_argument('--gap-penalty',
		dest = 'gap_penalty',
		default = 2,
		metavar = '',
		type = int,
		help = "Gap penalty"
	)
	issr_group.add_argument('--min-score',
		dest = 'min_required_score',
		default = 10,
		metavar = '',
		type = int,
		help = "minimum required score"
	)

	vntr_group = parser_search.add_argument_group('VNTR', 'VNTR search parameter')
	vntr_group.add_argument('--min-motif-len',
		dest = 'min_motif_length',
		default = 7,
		metavar = '',
		type = int,
		help = "minimum motif length"
	)
	vntr_group.add_argument('--max-motif-len',
		dest = 'max_motif_length',
		default = 30,
		metavar = '',
		type = int,
		help = "maximum motif length"
	)
	vntr_group.add_argument('--min-rep',
		dest = 'min_repeats',
		default = 2,
		metavar = '',
		type = int,
		help = "minimum repeats"
	)

	#extract flanking sequence parameters
	parser_flank = subparsers.add_parser('flank',
		help = "extract flanking sequence",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_flank.add_argument('-i', '--in',
		dest = 'infile',
		metavar = 'tsv',
		required = True,
		help = "output file of search command"
	)
	parser_flank.add_argument('-f', '--fasta',
		dest = 'seqfile',
		metavar = 'fasta',
		required = True,
		help = 'fasta formatted file containing sequences'
	)
	parser_flank.add_argument('-l', '--length',
		dest = 'length',
		metavar = 'int',
		required = True,
		default = 100,
		help = 'flanking sequence length'
	)
	parser_flank.add_argument('-o', '--out',
		dest = 'outfile',
		metavar = 'file',
		default = 'stdout',
		help = "output file name"
	)

	parser_locate = subparsers.add_parser('mapping', help="mapping tandem repeats to gene (CDS, UTR etc.)")
	parser_locate.add_argument('-i', '--in',
		dest = 'infile',
		metavar = 'tsv',
		required = True,
		help = 'output file of search command'
	)
	parser_locate.add_argument('-a', '--annot',
		dest = 'annotfile',
		metavar = 'gtf',
		required = True,
		help = 'annotation file with gff or gtf format'
	)
	parser_locate.add_argument('-o', '--out',
		dest = 'outfile',
		metavar = 'file',
		default = 'stdout',
		help = 'output file name (default: stdout)'
	)

	parser_primer = subparsers.add_parser('primer',
		help='Design primer for tandem repeats',
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_primer.add_argument('-i', '--in',
		dest = 'infile',
		metavar = 'tsv',
		required = True,
		help = "output file of search command"
	)
	parser_primer.add_argument('-f', '--fasta',
		dest = 'seqfile',
		metavar = 'fasta',
		required = True,
		help = 'fasta formatted file containing sequences'
	)
	parser_primer.add_argument('-o', '--out',
		dest = 'outfile',
		metavar = 'file',
		default = 'stdout',
		help = 'output file name (default: stdout)'
	)
	parser_primer.add_argument('-s', '--setting',
		dest = 'primer3_setting_file',
		required = True,
		metavar = 'setting file',
		help = 'primer3 setting file'
	)

	primer3_group = parser_primer.add_argument_group('Primer3', 'primer3 parameter')
	primer3_group.add_argument('--setting-file',
		dest = 'setting_file',
		required = True,
		metavar = 'file',
		help = 'primer3 setting file'
	)

	args = parser.parse_args()

	if args.cmd == 'search':
		if args.ssr_type == 'ssr':
			SSRSearchJob(args)
		
		elif args.ssr_type == 'cssr':
			CSSRSearchJob(args)
