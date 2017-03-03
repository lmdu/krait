#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import math
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from config import CACHE_PATH

def cache_path(name):
	return os.path.join(CACHE_PATH, "%s.png" % name)

class Plots:
	@staticmethod
	def SSRLenPie(data):
		x = []
		labels = []
		for row in data[1:]:
			x.append(row[1])
			labels.append(row[0])
		plt.pie(x, labels = labels, autopct = '%1.1f%%')
		plt.savefig(cache_path('ssr_len_pie_plot'), dpi=300, bbox_inches='tight')
		plt.close()

	@staticmethod
	def mostMoitfBar(data):
		motifs = []
		counts = []
		ra = []
		colors = []
		motif_len = 1
		temp_rows = []
		for row in sorted(data[1:], key=lambda x: (len(x[0]), -x[1])):
			if motif_len == len(row[0]):
				temp_rows.append(row)

			else:
				motif_len = len(row[0])
				if len(temp_rows) > 5:
					temp_rows = temp_rows[:5]
				
				for temp_row in temp_rows:
					motifs.append(temp_row[0])
					counts.append(temp_row[1])
					ra.append(temp_row[4])
					colors.append(len(row[0]))

				temp_rows = [row]

		if len(temp_rows) > 5:
			temp_rows = temp_rows[:5]
		for temp_row in temp_rows:
			motifs.append(temp_row[0])
			counts.append(temp_row[1])
			ra.append(temp_row[4])
			colors.append(len(row[0]))
		
		ind = range(len(motifs))
		plt.bar(ind, counts)
		plt.ylabel('SSR counts')
		plt.xlabel('Motifs')
		plt.xticks(ind, motifs, rotation=60, size='x-small')
		plt.savefig(cache_path('ssr_abundant_motif_bar_plot'), dpi=300, bbox_inches='tight')
		plt.close()

	@staticmethod
	def SSRMotifScatter(data):
		x = []
		y = []
		s = []
		c = []
		for row in data[1:]:
			x.append(math.log10(row[5]))
			y.append(math.log10(row[6]))
			s.append(row[4])
			c.append(len(row[0]))
		
		plt.scatter(x, y, s, c, alpha=0.6)

		
		plt.xlabel("log relative abundance")
		plt.ylabel("log relative density")
		#plt.legend()

		plt.savefig(cache_path('ssr_motif_scatter_plot'), dpi=300, bbox_inches='tight')
		plt.close()



		