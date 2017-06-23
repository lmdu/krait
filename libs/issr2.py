#!/usr/bin/env python

def print_matrix(m):
	for i in m:
		print " ".join(map(lambda x: "%-2d" % x,i))
	print "------------------------------------------"

def direction(left, top, lefttop):
	m = min((left, top, lefttop))
	if m == lefttop:
		return 0
	elif m == top:
		return 1
	else:
		return 2

def edit_distance(seq, motif, start):
	m = len(motif)
	N = 15
	size = len(seq) - start

	if size > N:
		size = N

	d = [[-1]*N for i in range(N)]

	for i in range(N):
		d[i][0] = i

	for j in range(N):
		d[0][j] = j

	error = 0
	matches = 0
	substitution = 0
	deletion = 0
	insertion = 0

	i = 1 #column no. in row
	j = 1 #row no. in column
	x = 1 #row position
	y = 1 #col position
	while(y<=size):
		#fill row, column number fixed
		if i != y:
			base = seq[start + y -1]
			for i in range(1, x):
				if base == motif[(i-1)%m]:
					d[i][y] = d[i-1][y-1]
				else:
					d[i][y] = min(d[i-1][y-1], d[i-1][y], d[i][y-1]) + 1

		#fill column, row number fixed
		if j != x:
			base = motif[(x-1)%m]
			for j in range(1, y):
				if base == seq[start + j -1]:
					d[x][j] = d[x-1][j-1]
				else:
					d[x][j] = min(d[x-1][j-1], d[x-1][j], d[x][j-1]) + 1
		
		i = y
		j = x

		if seq[start+y-1] == motif[(x-1)%m]:
			d[x][y] = d[x-1][y-1]
		else:
			d[x][y] = min(d[x-1][y-1], d[x-1][y], d[x][y-1]) + 1
		
			smaller = min(d[x][y], d[x-1][y], d[x][y-1])
			if smaller == d[x][y]:
				pass
			elif smaller == d[x][y-1]:
				y -= 1
			else:
				x -= 1

		x += 1
		y += 1
	
	#print matches, substitution, deletion, insertion
	print_matrix(d)
	
	i = x - error
	j = y - error
	while i>0 or j>0:
		direct direction(d[i-1][j], d[])
		if :
			i -= 1
			j -= 1


	print matches, substitution, insertion, deletion

if __name__ == '__main__':
	edit_distance('AAGAGAAG', 'AAG', 0)