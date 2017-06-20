#!/usr/bin/env python

def print_matrix(m):
	for i in m:
		print " ".join(map(lambda x: "%-2d" % x,i))

def edit_distance(seq, motif, start):
	m = len(motif)
	N = 200

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
	for k in range(1, N):
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
	while i>0 and j>0:
		if d[i][j] == d[i-1][j-1] and d[i-1][j-1]<=d[i-1][j] and d[i-1][j-1]<=d[i][j-1]:
			matches += 1
			i -= 1
			j -= 1
			continue

		prev = d[i][j] - 1

		if prev == d[i-1][j-1]:
			substitution += 1
			i -= 1
			j -= 1
		elif prev == d[i][j-1]:
			insertion += 1
			j -= 1
		else:
			deletion += 1
			i -= 1

	print matches, substitution, insertion, deletion


def	edDistRecursive(x,	y):
	if	len(x)	==	0:	return	len(y)
	if	len(y)	==	0:	return	len(x)
	delt	=	1	if	x[-1]	!=	y[-1]	else	0
	diag	=	edDistRecursive(x[:-1],	y[:-1])	+	delt
	vert	=	edDistRecursive(x[:-1],	y)	+	1
	horz	=	edDistRecursive(x,	y[:-1])	+	1
	return	min(diag,	vert,	horz)

if __name__ == '__main__':
	edit_distance('AAGAAGAAGACAGAAGAGAAGAAGACGAAGTTTT', 'AAG', 0)