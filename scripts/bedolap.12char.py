import sys
import os 
import yaml

with open("/diskmnt/Projects/ICGC_MC3/mc3_icgc_variant_pipeline/config.yaml", 'r') as ymlfile:
    config = yaml.load(ymlfile)


egmap = sys.argv[1]
ename = sys.argv[2]
gname = sys.argv[3]
enf = 0
gnf = 0


EXOME = {}
GENOM = {}


GENOME_ALIQUOT_TO_PCAWG_DONORID = {}
with open(config['GENOME_ALIQUOT_INFO']) as f:
    header = next(f)
    for line in f:
        _, __, submitter_donor_id, donor_id, _____, aliquot_id, *_ = line.rstrip().split('\t')
        GENOME_ALIQUOT_TO_PCAWG_DONORID[aliquot_id] = donor_id


# Set up EXOME MAPPING FOR WIGS 
EXOME_ALIQUOT_TO_ID = {}
with open(config['EXOME_ALIQUOT_INFO']) as f:
    header = next(f)
    for line in f:
        _, __, aliquot_id, ____, _____, exome_id, *_ = line.rstrip().split('\t')
        EXOME_ALIQUOT_TO_ID[aliquot_id] = exome_id

# Set up the EXOME 2 GENOME MAPPING
EXOME_SAMPLE_IDS_TO_GENOME = {}
GENOME_SAMPLE_IDS_TO_EXOME = {}

with open(config['OVERLAPPED_SAMPLE_INFO']) as f:
    header = next(f)
    for line in f:
        pcawg_aliquot, mc3_aliquot, exome_id, genome_id = line.rstrip().split('\t')
 #       if exome_id in SUBSETS: 
        EXOME_SAMPLE_IDS_TO_GENOME[exome_id] = GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]
        GENOME_SAMPLE_IDS_TO_EXOME[GENOME_ALIQUOT_TO_PCAWG_DONORID[pcawg_aliquot]] = exome_id



#with open(egmap,"r") as f:
#	eghead = f.readline().strip()
#	for line in f:
#		exome = line.strip().split(" ")[0]
#		genom = line.strip().split(" ")[1]
#		if exome[:12] == genom[:12]:
#			EXOME[exome] = genom
#			GENOM[genom] = exome
#f.close()


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
		if etcga in EXOME_SAMPLE_IDS_TO_GENOME:
			k = '_'.join([str(x) for x in [chrom,start,stop,etcga]])
			gtcga = EXOME_SAMPLE_IDS_TO_GENOME.get(etcga)
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
		gtcga = l[42]
		if gtcga in GENOME_SAMPLE_IDS_TO_EXOME:
			k = '_'.join([str(x) for x in [chrom,start,stop,gtcga]])
			etcga = GENOME_SAMPLE_IDS_TO_EXOME.get(gtcga)
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




