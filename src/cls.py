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
	def __init__(self, target, items, cpus):
		self.target = target
		self.items = items
		self.results = {}
		
		#CPU counts specified by user
		self.cpus = cpus
		
		#number of current tasks in pool
		self.tasks = 0

		if cpus > mp.cpu_count():
			self.cpus = mp.cpu_count()

		self.pool = mp.Pool(self.cpus)
		self.run_task()

	def success(self, res):
		self.tasks -= 1
		if res:
			self.results[res[0][0]] = res

	def failure(self, error):
		print("Error occured: {}".format(error))

	def add_task(self, args):
		self.pool.apply_async(
			func = self.target,
			args = args, 
			callback = self.success,
			error_callback = self.failure
		)
		self.tasks += 1

	def run_task(self):
		while 1:
			if self.tasks < self.cpus:
				try:
					item = next(self.items)
					self.add_task(item)
				
				except StopIteration:
					break
			else:
				time.sleep(0.001)
		
		self.pool.close()
		self.pool.join()

	def merge_task(self):
		with open(self.out, 'w') as fw:
			for outfile in glob.glob('{}.*'.format(self.out)):
				with open(outfile) as fh:
					fw.write(fh.read())
					os.remove(outfile)

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
	if outfile == 'stdout':
		fw = sys.stdout
	else:
		fw = open(outfile, 'w', newline='')

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

def search_cssr():
	pass

def search_issr():
	pass

def search_vntr():
	pass

def search_parameters(seqs, args):
	for name, seq in seqs:
		if args.ssr_type == 'ssr':
			yield (name, seq, args.repeats, args.level)
		
		elif args.ssr_type == 'cssr':
			yield (name, seq, args.repeats, args.dmax)

		elif args.ssr_type == 'issr':
			pass

		elif args.ssr_type == 'vntr':
			pass

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

	description = '''
	Search for microsatellite from large genome

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
		version = '%(prog)s version 1.0.0'
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
	parser_search.set_defaults(func=search_tandem)

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
		help = "input fasta file or gzipped fasta file"
	)
	parser_search.add_argument('-o', '--out',
		dest = 'outfile',
		default = 'stdout',
		metavar = '',
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

	parser_locate = subparsers.add_parser('mapping', help="mapping microsatellites or tandem repeats to gene (CDS, UTR etc.)")
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

	args.func(args)
