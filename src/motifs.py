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

def cycle_equal_motifs(motif):
	'''
	remove the similarity motifs e.g. AC and CA, CA is AC.
	@motif string.
	@return string, return the first motif.
	'''
	if len(motif) == 1:
		return [motif]
	
	return ["{}{}".format(motif[i+1:], motif[0:i+1]) for i, b in enumerate(motif)]

def complement_motifs(motif):
	'''
	motif completement eg. TG -> AC
	@para motif str
	@return list
	'''
	codes = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	complement = "".join([codes.get(b, b) for b in motif])
	return cycle_equal_motifs(complement)

def reverse_motifs(motif):
	return cycle_equal_motifs(motif[::-1])

def reverse_complement_motifs(motif):
	'''
	remove the reverse complete motifs e.g. AC and GT.
	@motif string
	@return list
	'''
	return complement_motifs(motif[::-1])

def motif_to_number(motif):
	sort_rule = {'A': '1', 'T': '2', 'C': '3', 'G': '4'}
	return int("".join(sort_rule.get(a, '5') for a in motif.upper()))

def motif_sorted(motifs):
	return sorted(motifs, key=motif_to_number)


class MotifStandard:
	_motif_mapping = {}

	def __init__(self, level=0):
		self.set_level(level)

	def set_level(self, level=0):
		'''Set standard level and generate motif mapping'''
		self.level = level
		if self.level > 0:
			self.generate_mapping()

	def get_standard(self, motif):
		'''Get standard motif for given motif'''
		if self.level == 0:
			return motif
		
		return self._motif_mapping.get(motif, None)

	def generate_mapping(self):
		'''Generate motif reverse mapping'''
		motifs = self.generate_motifs()
		self._motif_mapping = {m2:m1 for m1 in motifs for m2 in motifs[m1]}

	def generate_motifs(self):
		'''Generate standard motifs and corresponding equal motifs'''
		bases = ['A', 'T', 'C', 'G']
		motifs = []
		for i in range(1, 7):
			for motif in itertools.product(bases, repeat=i):
				motif = "".join(motif)
				if is_motif(motif):	
					motifs.append(motif)

		mappings = {}

		if self.level == 0:
			#5356 motifs
			return motifs

		for motif in motifs:
			candiates = []

			if self.level >= 1:
				candiates.extend(cycle_equal_motifs(motif))

			if self.level >= 2:
				candiates.extend(reverse_complement_motifs(motif))

			if self.level >= 3:
				candiates.extend(complement_motifs(motif))

			if self.level >= 4:
				candiates.extend(reverse_motifs(motif))

			#remove the same motifs
			candiates = list(set(candiates))

			#sort motifs as A>T>C>G
			candiates = motif_sorted(candiates)

			if candiates[0] not in mappings:
				mappings[candiates[0]] = candiates

		return mappings

if __name__ == '__main__':
	s = MotifStandard(2)
	print(s.get_standard('ATG'))
