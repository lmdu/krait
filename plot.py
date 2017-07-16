#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')

from config import CACHE_PATH

def cache_path(name):
	return os.path.join(CACHE_PATH, "%s.png" % name)

def line(x, y, name, dpi=96):
	pass

def pie(x, label, name, dpi=96):
	plt.pie(x, labels=label, autopct='%1.1f%%')
	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()

def bar(x, y, xlabel, ylabel, name, dpi=96):
	ind = range(len(x))
	plt.bar(ind, y)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.xticks(ind, x, rotation=60, size='x-small')
	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()

def box(x, labels, name, dpi=96):
	colors = ['pink', 'lightblue', 'lightgreen', 'red', 'yellow', 'purple']
	bplot = plt.boxplot(x, 
		notch = True, 
		vert = True, 
		patch_artist = True,
		labels = labels,
		showfliers = False
	)
	#for path, color in zip(bplot['boxes'], colors):
	#	path.set_facecolor(color)

	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()


