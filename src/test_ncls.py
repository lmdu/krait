from ncls import NCLS
import numpy as np

starts = []
ends = []
ids = []
with open('data.txt') as fh:
	for line in fh:
		cols = line.strip().split()
		starts.append(int(cols[0]))
		ends.append(int(cols[1]))
		ids.append(int(cols[2]))

		if int(cols[2]) > 70000:
			break

starts = np.array(starts, dtype=np.long)
ends = np.array(ends, dtype=np.long)
ids = np.array(ids, dtype=np.long)

ncls = NCLS(starts, ends, ids)

for i in ncls.find_overlap(76623690, 76624000):
	print(i)
