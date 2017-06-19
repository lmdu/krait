#!/usr/bin/env python

def print_matrix(m):
	for i in m:
		print " ".join(map(str,i))

def edit_distance(seq, motif, start):
	m = len(motif)
	N = 10

	d = [[0]*N for i in range(N)]

	for i in range(N):
		d[i][0] = i

	for j in range(N):
		d[0][j] = j

	for i in range(1, N):
		for j in range(1, N):
			if seq[start+i] == motif[j%m]:
				d[i][j] = d[i-1][j-1]
			else:
				d[i][j] = min(d[i-1][j-1], d[i-1][j], d[i][j-1]) + 1

			


	

if __name__ == '__main__':
	edit_distance('abc', 'ab', 0)