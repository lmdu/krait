import sys
import argparse
import threading
import multiprocessing

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
	required = True,
	help = "input fasta file or gzipped fasta file"
)
parser_search.add_argument('-o', '--out',
	dest = 'outfile',
	default = 'stdout',
	metavar = '',
	help = "output file name"
)
parser_search.add_argument('-f', '--outfmt',
	dest = 'outfmt',
	default = 'tsv',
	choices = ('tsv', 'csv', 'gff'),
	metavar = '',
	help="output format [csv, tsv, gff]"
)
parser_search.add_argument('-c', '--cpus',
	dest = 'cpus',
	default = 1,
	type = int,
	metavar = '',
	help = 'number of threads'
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

class ProcessPool:
	def __init__(self):
		pass
