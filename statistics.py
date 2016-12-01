def getTotalBases(pf):
	return sum([ len(pf[name]) for name in pf ])

def getSeqCount(pf):
	return len(pf.keys())

class Statistics:
	def __init__(self, db, work="ssrhit"):
		self.db = db
		self.work = work
		self.groups = {
			1: 'mononucleotide', 
			2: 'dinucleotide', 
			3: 'trinucleotide', 
			4: 'tetranucleotide',
			5: 'pentanucleotide',
			6: 'hexanucleotide'
		}

	def do(self):
		self.summaryStat()
		self.motifStat()
		self.repeatStat()
		self.seqStat()

	def summaryStat(self):
		fp = open("%s_summary_stat.xls" % self.work, "w")
		totalSSRs = self.db.get("SELECT COUNT(1) FROM ssr")
		totalSSRLen = self.db.get("SELECT SUM(length) FROM ssr")
		totalSeqs = self.db.getStat("seqNum")
		totalSeqLen = self.db.getStat("baseNum")
		avgLen = round(totalSSRLen/totalSSRs, 2)
		frequency = round(totalSSRs/(totalSeqLen/1000000), 2)
		density = round(totalSSRLen/(totalSeqLen/1000000), 2)
		fp.write("Total Sequences:\t%s\n" % totalSeqs)
		fp.write("Total Sequence Length:\t%s\n" % totalSeqLen)
		fp.write("Total SSRs:\t%s\n" % totalSSRs)
		fp.write("SSR Average Length:\t%s\n" % avgLen)
		fp.write("SSR Frequency (loci/Mb):\t%s\n" % frequency)
		fp.write("SSR Density (bp/Mb):\t%s\n" % density)

		sql = (
			"SELECT length(motif) AS m, COUNT(1) AS c "
			"SUM(length) AS l FROM ssr GROUP BY m"
		)
		
		fp.write("\nClass\tCount\tLength\n")
		
		for row in self.db.exe(sql):
			fp.write("%s\t%s\n" % (self.groups[row[0]], row[1], row[2]))

		fp.close()

	def motifStat(self):
		sql = (
			"SELECT motif2, COUNT(1), SUM(length) FROM ssr "
			"GROUP BY motif2 ORDER BY length(motif2),motif2"
		)
		with open("%s_motif_stat.xls" % self.work, "w") as fp:
			fp.write("Motif\tCount\tLength\n")
			for row in self.db.exe(sql):
				fp.write("%s\t%s\n" % row)

	def repeatStat(self):
		sql = (
			"SELECT length(motif) AS m, repeat, COUNT(1) AS c "
			"FROM ssr GROUP BY m, repeat"
		)
		with open("%s_repeat_stat.xls" % self.work, "w") as fp:
			fp.write("Class\tRepeat\tCount\n")
			for row in self.db.exe(sql):
				fp.write("%s\t%s\t%s\n" % (self.groups[row[0]], row[1], row[2]))

	def seqStat(self):
		sql = (
			"SELECT chrom, COUNT(1) AS c, SUM(length) AS srrl, "
			"value AS seql FROM ssr,stat WHERE chrom=name "
			"GROUP BY chrom ORDER BY c DESC"
		)
		with open("%s_seq_stat.xls" % self.work, "w") as fp:
			fp.write("Seq Name\tSeq Length\tSSR Count\tSSR Length\n")
			for row in self.db.exe(sql):
				fp.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % row)