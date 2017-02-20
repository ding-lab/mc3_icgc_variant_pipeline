#!/gscuser/sfoltz/anaconda3/bin/python

import sys
import os 

egmap = sys.argv[1]
ename = sys.argv[2]
gname = sys.argv[3]
enf = 0
gnf = 0


EXOME = {}
GENOM = {}

with open(egmap,"r") as f:
	eghead = f.readline().strip()
	for line in f:
		exome = line.strip().split(" ")[0]
		genom = line.strip().split(" ")[1]
		if exome[:12] == genom[:12]:
			EXOME[exome] = genom
			GENOM[genom] = exome
f.close()


EMAF = {}
GMAF = {}

with open(ename,"r") as f:
	for line in f:
		l = line.strip().split("\t")
		enf = len(l)
		chrom = l[0]
		start = l[1]
		stop = l[2]
		etcga = l[11]
		if etcga in EXOME:
			k = '_'.join([str(x) for x in [chrom,start,stop,etcga]])
			gtcga = EXOME.get(etcga)
			v = '_'.join([str(x) for x in [chrom,start,stop,gtcga]])
			elist = [v,line.strip()]
			EMAF[k] = elist
f.close()


with open(gname,"r") as f:
	for line in f:
		l = line.strip().split("\t")
		gnf = len(l)
		chrom = l[0]
		start = l[1]
		stop = l[2]
		gtcga = l[45]
		if gtcga in GENOM:
			k = '_'.join([str(x) for x in [chrom,start,stop,gtcga]])
			etcga = GENOM.get(gtcga)
			v = '_'.join([str(x) for x in [chrom,start,stop,etcga]])
			elist = [v,line.strip()]
			GMAF[k] = elist
f.close()


emissing = "\t".join("NA" for i in range(enf))
gmissing = "\t".join("NA" for i in range(gnf))


for item in EMAF:
	newkey = EMAF[item][0]
	eline = EMAF[item][1]
	if newkey in GMAF:
		print(eline+"\t"+str(GMAF[newkey][1]))
		del GMAF[newkey]
	else: 
		print(eline+"\t"+gmissing)

for item in GMAF:
	gline = GMAF[item][1]
	print(emissing+"\t"+gline)




