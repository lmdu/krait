#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import math
import matplotlib
from numpy import sin, cos, pi
matplotlib.use('Agg')
matplotlib.rc('font', size=8)
import matplotlib.pyplot as plt
plt.style.use('seaborn-paper')

from config import CACHE_PATH

def cache_path(name):
	return os.path.join(CACHE_PATH, "%s.png" % name)

def line(x, y, xlabel, ylabel, name, dpi=96):
	plt.plot(x, y, linewidth=2)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()

def pie(x, label, name, dpi=96):
	data = sorted(zip(x, label), key=lambda x:x[0], reverse=True)
	x = []
	label = []
	for a, b in data:
		x.append(a)
		label.append(b)

	total = sum(x)
	fracs = [round(v/total*100, 2) for v in x]
	p = plt.pie(fracs, startangle=90)
	plt.axis('equal')
	#Annotations from https://commons.wikimedia.org/wiki/File:PV-Global-2010.svg
	y1s = []
	for p1, l1, f1 in zip(p[0], label, fracs):
		r = p1.r
		dr = r*0.1
		t1, t2 = p1.theta1, p1.theta2
		theta = (t1+t2)/2.

		xc, yc = r/2.*cos(theta/180.*pi), r/2.*sin(theta/180.*pi)
		x1, y1 = (r+dr)*cos(theta/180.*pi), (r+dr)*sin(theta/180.*pi)
		if x1 > 0 :
			x1 = r+2*dr
			ha, va = "left", "center"
			tt = -180
			cstyle = "angle,angleA=0,angleB=%f" % (theta,)
		else:
			x1 = -(r+2*dr)
			ha, va = "right", "center"
			tt = 0
			cstyle = "angle,angleA=0,angleB=%f" % (theta,)

		# Prevent overlapping text
		y1 = round(y1*10)/10
		if any([(abs(y1 - y) < 0.15) for y in y1s]):
			y1 = max(y1s) + 0.15

		l1 = l1 + ': ' + str(f1) + '%'
		plt.annotate(l1, (xc, yc), 
			xycoords = "data",
			xytext = (x1, y1),
			textcoords = "data",
			ha = ha,
			va = va,
			arrowprops=dict(arrowstyle="-", connectionstyle=cstyle, patchB=p1)
		)

		y1s.append(y1)


	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()

def bar(x, y, xlabel, ylabel, name, dpi=96):
	n = range(len(x))
	plt.bar(n, y)
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.xticks(n, x, rotation=90, size='x-small')
	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()

def box(x, xlabels, ylabel, name, dpi=96):
	colors = ['pink', 'lightblue', 'lightgreen', 'orangered', 'goldenrod', 'orchid']
	bplot = plt.boxplot(x,
		notch = False,
		vert = True, 
		patch_artist = True,
		labels = xlabels,
		showfliers = False
	)
	for path, color in zip(bplot['boxes'], colors):
		path.set_facecolor(color)
	plt.ylabel(ylabel)
	plt.savefig(cache_path(name), dpi=dpi, bbox_inches='tight')
	plt.close()


