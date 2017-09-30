#!/usr/bin/env python
import pyfaidx

def initial_matrix(size):
	matrix = []
	for i in range(size+1):
		matrix.append([])
		for j in range(size+1):
			if i == 0:
				matrix[i].append(j)
			elif j == 0:
				matrix[i].append(i)
			else:
				matrix[i].append(0)

	return matrix

def print_matrix(matrix):
	for i in range(len(matrix)):
		print "\t".join(map(str,matrix[i]))
		
def build_left_matrix(seq, motif, matrix, start, size, max_error):
	ref1 = None
	ref2 = None
	i = 1
	j = 1
	x = 1
	y = 1
	last_x = 0
	last_y = 0
	mlen = len(motif)
	error = 0
	smaller = 0
	
	while y<=size:
		ref1 = seq[start-y]
		ref2 = motif[(mlen-x%mlen)%mlen]
		#fill column column number fixed
		if i != y:
			for i in range(1, x):
				if ref1 == motif[(mlen-i%mlen)%mlen]:
					matrix[i][y] = matrix[i-1][y-1]
				else:
					matrix[i][y] = min(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1

		#fill row, row number fixed
		if j != x:
			for j in range(1, y):
				if ref2 == seq[start-j]:
					matrix[x][j] = matrix[x-1][j-1];
				else:
					matrix[x][j] = min(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1

		i = y
		j = x

		if ref1 == ref2:
			matrix[x][y] = matrix[x-1][y-1]
			error = 0
		else:
			if error == 0:
				last_x = x - 1
				last_y = y - 1
			
			error += 1

			if error > max_error:
				break
			
			matrix[x][y] = min(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1
		
		smaller = min(matrix[x][y], matrix[x-1][y], matrix[x][y-1])

		if smaller != matrix[x][y]:
			if matrix[x-1][y] != matrix[x][y-1]:
				if smaller == matrix[x][y-1]:
					y -= 1
				else:
					x -= 1

		x += 1
		y += 1

	#print x, y

	#print_matrix(matrix)
			
	if error:
		res = (last_x, last_y)
	else:
		res = (x-1, y-1)

	return res

def build_right_matrix(seq, motif, matrix, start, size, max_error):
	ref1 = None
	ref2 = None
	i = 1
	j = 1
	x = 1
	y = 1
	last_x = 0
	last_y = 0
	mlen = len(motif)
	error = 0
	smaller = 0
	
	while y<=size:
		ref1 = seq[start+y]
		ref2 = motif[(x-1)%mlen]
		#fill column column number fixed
		if i != y:
			for i in range(1, x):
				if ref1 == motif[(i-1)%mlen]:
					matrix[i][y] = matrix[i-1][y-1]
				else:
					matrix[i][y] = min(matrix[i-1][y-1], matrix[i-1][y], matrix[i][y-1]) + 1

		#fill row, row number fixed
		if j != x:
			for j in range(1, y):
				if ref2 == seq[start+j]:
					matrix[x][j] = matrix[x-1][j-1];
				else:
					matrix[x][j] = min(matrix[x-1][j-1], matrix[x-1][j], matrix[x][j-1]) + 1

		i = y
		j = x

		if ref1 == ref2:
			matrix[x][y] = matrix[x-1][y-1]
			error = 0
		else:
			if error == 0:
				last_x = x-1
				last_y = y-1
			
			error += 1

			if error > max_error:
				break
			
			matrix[x][y] = min(matrix[x-1][y-1], matrix[x-1][y], matrix[x][y-1]) + 1
		
		smaller = min(matrix[x][y], matrix[x-1][y], matrix[x][y-1])

		if smaller != matrix[x][y]:
			if matrix[x-1][y] != matrix[x][y-1]:
				if smaller == matrix[x][y-1]:
					y -= 1
				else:
					x -= 1

		x += 1
		y += 1
			
	if error:
		res = (last_x, last_y)
	else:
		res = (x-1, y-1)

	return res


def backtrace_matrix(matrix, diagonal):
	i, j = diagonal

	substitution = 0
	insertion = 0
	deletion = 0
	match = 0
	cost = 0
	r = j

	while i>0 and j>0:
		cost = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1])
		if cost == matrix[i-1][j-1]:
			if cost == matrix[i][j]:
				match += 1
			else:
				substitution += 1

			i -= 1
			j -= 1

		elif cost == matrix[i-1][j]:
			deletion += 1
			i -= 1

		else:
			insertion += 1
			j -= 1

	if i>0:
		deletion += 1
	elif j>0:
		insertion += 1

	return (r, match, substitution, insertion, deletion)


def search_issr(seq, seed_repeats, seed_minlen, max_errors, mis_penalty, gap_penalty, required_score, size):
	res = []
	matrix = initial_matrix(size)
	seqlen = len(seq)
	i = 0
	while i < seqlen:
		if seq[i] == 'N':
			i += 1
			continue

		j = 1
		while j <= 6:
			seed_start = i
			seed_length = j
			while seed_start+seed_length < seqlen and seq[i] == seq[i+j] and seq[i] != 'N':
				i += 1
				seed_length += 1

			seed_repeat = seed_length/j

			if seed_repeat >= seed_repeats and seed_length >= seed_minlen:
				motif = seq[seed_start:seed_start+j]

				seed_end = seed_start + seed_length - seed_length%j - 1
				matches = seed_length - seed_length%j
				insertion = 0
				deletion = 0
				substitution = 0
				
				#extend left
				
				extend_start = seed_start
				extend_len = extend_start
				if extend_len > size:
					extend_len = size
				extend_end = build_left_matrix(seq, motif, matrix, extend_start, extend_len, max_errors)
				extend_ok = backtrace_matrix(matrix, extend_end)
				start = extend_start - extend_ok[0] + 1
				matches += extend_ok[1]
				substitution += extend_ok[2]
				insertion += extend_ok[3]
				deletion += extend_ok[4]
				
				#extend right
				extend_start = seed_end
				extend_len = seqlen - extend_start - 1
				if extend_len > size:
					extend_len = size
				extend_end = build_right_matrix(seq, motif, matrix, extend_start, extend_len, max_errors)
				extend_ok = backtrace_matrix(matrix, extend_end)
				end = extend_start + extend_ok[0] + 1
				length = end - start + 1
				matches += extend_ok[1]
				substitution += extend_ok[2]
				insertion += extend_ok[3]
				deletion += extend_ok[4]

				score = matches - substitution*mis_penalty - (insertion+deletion)*gap_penalty
				
				if score >= required_score:
					res = (motif, j, start, end, length, matches, substitution, insertion, deletion, score)
					print "%s\t%s" % (",".join(map(str, res)), seq[start-1:end])
					i = end
					j = 0
				else:
					i = seed_start
			else:
				i = seed_start
			j += 1
		i += 1

	#return res



def left_alignment(seq, motif, start, matrix, diagonal):
	i, j = diagonal
	mlen = len(motif)
	origin = []
	copy = []
	while i>0 and j>0:
		ref1 = seq[start-j]
		ref2 = motif[(mlen-i%mlen)%mlen]
		cost = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1])
		if cost == matrix[i-1][j-1]:
			origin.append(ref1)
			copy.append(ref2)

			i -= 1
			j -= 1

		elif cost == matrix[i-1][j]:
			origin.append('-')
			copy.append(ref2)
			i -= 1

		else:
			origin.append(ref1)
			copy.append('-')
			j -= 1

	#print i, j

	if i>0:
		pass

	elif j>0:
		pass

	return (origin, copy)

def right_alignment(seq, motif, start, matrix, diagonal):
	i, j = diagonal
	mlen = len(motif)
	origin = []
	copy = []
	while i>0 and j>0:
		ref1 = seq[start+j]
		ref2 = motif[(i-1)%mlen]
		cost = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1])
		if cost == matrix[i-1][j-1]:
			origin.append(ref1)
			copy.append(ref2)

			i -= 1
			j -= 1

		elif cost == matrix[i-1][j]:
			origin.append('-')
			copy.append(ref2)
			i -= 1

		else:
			origin.append(ref1)
			copy.append('-')
			j -= 1

	if i>0:
		pass

	elif j>0:
		pass

	origin.reverse()
	copy.reverse()

	return (origin, copy)


def generate_alignment(seq, seed_repeats, seed_minlen, max_errors, size):
	matrix = initial_matrix(size)
	seqlen = len(seq)
	origin = []
	copy = []
	i = 0
	while i < seqlen:
		if seq[i] == 'N':
			i += 1
			continue

		j = 1
		while j <= 6:
			seed_start = i
			seed_length = j
			while seed_start+seed_length < seqlen and seq[i] == seq[i+j] and seq[i] != 'N':
				i += 1
				seed_length += 1

			seed_repeat = seed_length/j

			if seed_repeat >= seed_repeats and seed_length >= seed_minlen:
				motif = seq[seed_start:seed_start+j]

				seed_end = seed_start + seed_length - seed_length%j - 1
	
				#extend left	
				extend_start = seed_start
				extend_len = extend_start
				if extend_len > size:
					extend_len = size
				extend_end = build_left_matrix(seq, motif, matrix, extend_start, extend_len, max_errors)
				
				o, c = left_alignment(seq, motif, extend_start, matrix, extend_end)
				origin.extend(o)
				copy.extend(c)

				#seed alignment
				for k in range(seed_start, seed_end+1):
					origin.append(seq[k])
					copy.append(seq[k])
				
				#extend right
				extend_start = seed_end
				extend_len = seqlen - extend_start - 1
				if extend_len > size:
					extend_len = size
				extend_end = build_right_matrix(seq, motif, matrix, extend_start, extend_len, max_errors)
				o, c = right_alignment(seq, motif, extend_start, matrix, extend_end)
				origin.extend(o)
				copy.extend(c)
				
				return origin, copy
				
			else:
				i = seed_start
			j += 1
		i += 1


if __name__ == '__main__':
	#fasta = pyfaidx.Fasta('test.fa')
	#seq = str(fasta['NC_000913.3'])
	seq = "AAGATAAGAAGAAGATGAAGAGAAGTTTTTTT"
	#seq = "AAAAAAAAATCATTT"
	#seq = "TCATCATCAAAATCGCCAT"
	o, c = generate_alignment(seq, 3, 8, 3, 50)
	print "".join(o)
	print "".join(c)
	#seq = "GTTGTTGTTGATTG"
	#search_issr(seq, 3, 8, 3, 2, 5, 3, 20)

		
	
