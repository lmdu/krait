#!/usr/bin/env python
# -*- coding: utf-8 -*-
import itertools

def isMiniMotif(motif):
	'''
	check the motif length is weather or not minimal length,
	e.g. AAAA is mono- motif A not tetra- motif.
	@motif string.
	@return bool.
	'''
	if len(set(motif)) == 1: return True
	
	if len(motif) % 2 != 0: return False
	
	if len(motif) == 4:
		return True if motif[0:2] == motif[2:4] else False

	if len(motif) == 6:
		if motif[0:3] == motif[3:6]:
			return True
		elif motif[0:2] == motif[2:4] and motif[2:4] == motif[4:6]:
			return True
		else:
			return False

def getSimilarMotif(motif):
	'''
	remove the similarity motifs e.g. AC and CA, CA is AC.
	@motif string.
	@return string, return the first motif.
	'''
	if len(motif) == 2: return motif[::-1]
	similar_motifs = []
	for i, b in enumerate(motif):
		new_motif = "%s%s" % (motif[i+1:], motif[0:i+1])
		similar_motifs.append(new_motif)
		if new_motif == motif: break
	similar_motifs.sort()
	return similar_motifs[0]

def reverseCompleteMotif(motif):
	'''
	remove the reverse complete motifs e.g. AC and GT.
	@motif string.
	@return string, the reversed and completed motif.
	'''
	codes = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	return "".join([ codes[c] for c in motif[::-1]])

def generateMotifTable(mlen):
	'''
	generate all motifs by motif length.
	@mlen integer, motif length 2~7.
	@return dict, the motifs table.
	'''
	bases = ['A', 'C', 'G', 'T']
	table = {}
	for motif in itertools.product(bases, repeat=mlen):
		motif = "".join(list(motif))
		
		if isMiniMotif(motif): continue
		
		new_motif = getSimilarMotif(motif)
		if new_motif in table:
			table[new_motif].append(motif)
			continue

		new_motif = reverseCompleteMotif(new_motif)
		if new_motif in table:
			table[new_motif].append(motif)
			continue

		new_motif = reverseCompleteMotif(motif)
		if new_motif in table:
			table[new_motif].append(motif)
			continue

		new_motif = getSimilarMotif(new_motif)
		if new_motif in table:
			table[new_motif].append(motif)
			continue

		table[motif] = []
	
	return table

MOTIFS = {'A': 'A', 'C': 'C', 'T':'A', 'G':'C'}

#generate the all motifs 1~6.
for i in range(2, 7):
	motifs = generateMotifTable(i)
	for m1 in motifs:
		MOTIFS[m1] = m1
		for m2 in motifs[m1]:
			MOTIFS[m2] = m1

def normalize(motif):
	return MOTIFS[motif]



if __name__ == '__main__':
	import json
	with open("motif.json", "w") as fp:
		json.dump(MOTIFS, fp, indent=4)