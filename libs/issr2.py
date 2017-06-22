#!/usr/bin/env python

def print_matrix(m):
	for i in m:
		print " ".join(map(lambda x: "%-2d" % x,i))

def edit_distance(seq, motif, start):
	m = len(motif)
	N = 20
	size = len(seq) - start

	if size > N:
		size = N

	d = [[0]*N for i in range(N)]

	for i in range(N):
		d[i][0] = i

	for j in range(N):
		d[0][j] = j

	error = 0
	matches = 0
	substitution = 0
	deletion = 0
	insertion = 0
	for k in range(1, size):
		for l in range(1, k):
			#fill column, k is row
			if seq[start+l-1] == motif[(k-1)%m]:
				d[k][l] = d[k-1][l-1]
			else:
				d[k][l] = min(d[k-1][l-1], d[k-1][l], d[k][l-1]) + 1

			#fill row, k is column
			if seq[start+k-1] == motif[(l-1)%m]:
				d[l][k] = d[l-1][k-1]
			else:
				d[l][k] = min(d[l-1][k-1], d[l-1][k], d[l][k-1]) + 1

		#fill diagonal
		if seq[start+k-1] == motif[(k-1)%m]:
			d[k][k] = d[k-1][k-1]
			error = 0
		else:
			d[k][k] = min(d[k-1][k-1], d[k-1][k], d[k][k-1]) + 1
			error += 1

		if error > 3:
			break

	print_matrix(d)

	m = k - error
	print m

	i = m
	j = m
	while i>0 or j>0:
		if d[i][j] == d[i][j-1]+1:
			j -= 1
			insertion += 1
		elif d[i][j] == d[i-1][j]+1:
			i -= 1
			deletion += 1
		else:
			if d[i][j] == d[i-1][j-1]:
				matches += 1
			else:
				substitution += 1

			i -= 1
			j -= 1

	print matches, substitution, insertion, deletion

if __name__ == '__main__':
	edit_distance('AAGAAGAGAAG', 'AAG', 0)