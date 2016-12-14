#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pyfaidx

class Data(dict):
	def __getattr__(self, name):
		try:
			return self[name]
		except KeyError:
			raise AttributeError(name)

	def __setattr__(self, name, val):
		self[name] = val


class SequenceHighlighter:
	html = '''
	<!DOCTYPE html>
	<html>
	<head>
		<meta charset="utf-8">
	</head>
	<style>
		*{
			font-family: "Microsoft YaHei";
			margin:0;
			padding:0;
		}
		body{
			width: 100%%;
			background: #eee;
		}
		span{
			display: inline-block;
			position: relative;
		}
		.meta{
			position: absolute;
			font-size: 12px;
			color: black;
			top: -17px;
			left: 0;
			min-width: 200px;
			width: 100%%;
			border: 1px solid red;
			text-align: left;
		}
		.A,.T,.C,.G,.N{
			font-size: 16px;
			font-weight: bold;
			color: #fff;
			width:18px;
			line-height: 24px;
			text-align: center;
			margin-top: 20px;
		}
		.A{
			border-color: #088ef0;
			background: -webkit-gradient(linear, left top, left bottom, from(#34a5f8), to(#088ef0));
  			background: linear-gradient(#34a5f8, #088ef0);
		}
		.T{
			border-color: #fea502;
			background: -webkit-gradient(linear, left top, left bottom, from(#feb734), to(#fea502));
			background: linear-gradient(#feb734, #fea502);
		}
		.C{
			border-color: #ff2939;
			background: -webkit-gradient(linear, left top, left bottom, from(#ff5c69), to(#ff2939));
			background: linear-gradient(#ff5c69, #ff2939);
		}
		.G{
			border-color: #9ad824;
			background: -webkit-gradient(linear, left top, left bottom, from(#afe24d), to(#9ad824));
			background: linear-gradient(#afe24d, #9ad824);
		}
		.N{
			border-color: #665ce6;
  			background: -webkit-gradient(linear, left top, left bottom, from(#9088ec), to(#665ce6));
  			background: linear-gradient(#9088ec, #665ce6); 
		}
		.ssr{
			display:inline;
			border-bottom: 4px solid #000;
		}
	</style>
	<body>
		%s
	<body>
	</html>
	'''
	def __init__(self):
		self.dna = []

	def format_base(self, base):
		return '<span class="{b}">{b}</span>'.format(b=base)

	def format_meta(self, base, meta):
		return '<span class="{b}">{b}<div class="meta">{m}</div></span>'.format(b=base,m=meta)
	
	def format_flank(self, flank, meta=None):
		if meta:
			self.dna.append(self.format_meta(flank[0], meta))
			self.dna.append("".join([self.format_base(base) for base in flank[1:]]))
		else:
			self.dna.append("".join([self.format_base(base) for base in flank]))

	def format_ssr(self, ssr):
		ssr = "".join([self.format_base(base) for base in ssr])
		self.dna.append('<span class="ssr">%s</span>' % ssr)

	def render(self):
		return self.html % "".join(self.dna)


def get_ssr_sequence(seq_file, seq_name, start, stop, flank):
	'''
	Get the SSR sequence and flanking sequences
	@para seq_file, the file path of the fasta sequence
	@para seq_name, the name of the fasta sequence
	@para start, the start position of SSR
	@para stop, the stop position of SSR
	@para flank, the length of the flanking sequence
	@return ssr sequence with flanking sequences
	'''
	fasta = pyfaidx.Fasta(seq_file, sequence_always_upper=True)
	
	#get ssr sequence
	ssr = fasta[seq_name][start-1:stop].seq
	
	#get left flanking sequence
	left_flank_start = start - flank - 1
	if left_flank_start < 0:
		left_flank_start = 0
	left_flank = fasta[seq_name][left_flank_start:start]
	
	seq_len = len(fasta[seq_name][:])
	
	#get right flanking sequence
	right_flank_stop = stop + flank
	if right_flank_stop > seq_len:
		right_flank_stop = seq_len
	right_flank = fasta[seq_name][stop:right_flank_stop]

	highlighter = SequenceHighlighter()
	meta = '%s:%s-%s %s' % (seq_name, left_flank_start+1, start, len(left_flank))
	highlighter.format_flank(left_flank, meta)
	highlighter.format_ssr(ssr)
	highlighter.format_flank(right_flank)
	return highlighter.render()
