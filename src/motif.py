# -*- coding: utf-8 -*-
import itertools

def is_motif(motif):
	'''
	check the motif length is weather or not minimal length,
	e.g. AAAA is mono- motif A not tetra- motif
	@motif string
	@return bool
	'''
	if len(motif) == 1:
		return True

	if len(set(motif)) == 1:
		return False

	if len(motif) in (3,5):
		return True

	if len(motif) == 2:
		if motif[0] == motif[-1]:
			return False
		else:
			return True

	if len(motif) == 4:
		if motif[0:2] == motif[2:4]:
			return False
		else:
			return True

	if len(motif) == 6:
		if motif[0:3] == motif[3:6]:
			return False

		elif motif[0:2] == motif[2:4] == motif[4:6]:
			return False

		else:
			return True

def similar_motif(motif):
	'''
	remove the similarity motifs e.g. AC and CA, CA is AC.
	@motif string.
	@return string, return the first motif.
	'''
	if len(motif) == 1:
		return [motif]
	
	similar_motifs = []
	for i, b in enumerate(motif):
		new_motif = "%s%s" % (motif[i+1:], motif[0:i+1])
		similar_motifs.append(new_motif)
	return similar_motifs

def complete_motif(motif):
	'''
	motif completement eg. TG -> AG
	@para motif str
	@return list
	'''
	codes = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	new_motif = "".join([ codes.get(base, base) for base in motif])
	return similar_motif(new_motif)


def reverse_complete_motif(motif):
	'''
	remove the reverse complete motifs e.g. AC and GT.
	@motif string
	@return list
	'''
	return complete_motif(motif[::-1])

def reverse_motif(motif):
	return similar_motif(motif[::-1])

def motif_to_number(motif):
	sort_rule = {'A':'1', 'T':'2', 'C':'3', 'G':'4'}
	return int("".join(sort_rule.get(a, '5') for a in motif.upper()))

def motif_sorted(motifs):
	return sorted(motifs, key=motif_to_number)


class StandardMotif:
	_motifs = {}
	def __init__(self, level=0):
		self.level = level

	def setLevel(self, level=0):
		self._motifs = {}
		self.level = level

	def standard(self, motif):
		if self.level == 0:
			return motif

		if motif in self._motifs:
			return self._motifs[motif]

		motifs = []
		
		if self.level >= 1:
			motifs.extend(similar_motif(motif))

		if self.level >= 2:
			motifs.extend(reverse_complete_motif(motif))

		if self.level >= 3:
			motifs.extend(complete_motif(motif))

		if self.level >= 4:
			motifs.extend(reverse_motif(motif))


		#remove the same motifs
		motifs = list(set(motifs))

		#sort motifs as A>T>C>G
		motifs = motif_sorted(motifs)

		for motif in motifs:
			self._motifs[motif] = motifs[0]

		return motifs[0]

	def mapping(self):
		bases = ['A', 'T', 'C', 'G']
		motifs = {}
		for i in range(6):
			for motif in itertools.product(bases, repeat=i+1):
				motif = "".join(list(motif))
				if not is_motif(motif):
					continue

				smotif = self.standard(motif)
				if smotif not in motifs:
					motifs[smotif] = []

				if motif not in motifs[smotif]:
					motifs[smotif].append(motif)

		return motifs

