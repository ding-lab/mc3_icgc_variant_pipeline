import os
import sys

egmap = sys.argv[1]
ename = sys.argv[2]
gname = sys.argv[3]

EXOME = {}
GENOM = {}

with open(egmap,"r") as f:
    eghead = f.readline().strip()
    for line in f:
        exome = line.strip().split(" ")[0]
        genom = line.strip().split(" ")[1]
        EXOME[exome] = 1
        GENOM[genom] = 1
f.close()

with open("exome.maf", "w") as e:
    with open(ename, "r") as f:
        ehead = f.readline().strip().split("\t")
        hugo = ehead.pop(0)
        ehead.append(hugo)
        del ehead[:3]
        for item in ehead:
            e.write("%s\t" % item)
        e.write("\n")
        for line in f:
            eid = line.strip().split("\t")[15]
            l = line.strip().split("\t")
            gene = l.pop(0)
            l.append(gene)
            del l[:3]
            if eid in EXOME:
                for item in l:
                    e.write("%s\t" % item)
                e.write("\n")
    f.close()
e.close()


with open("genome.maf","w") as g:
    with open(gname,"r") as f:
        ghead = f.readline().strip().split("\t")
        hugo = ghead.pop(0)
        ghead.append(hugo)
        for item in ghead:
            g.write("%s\t" % item)
        g.write("\n")
        for line in f:
            gid = line.strip().split("\t")[46]
            l = line.strip().split("\t")
            gene = l.pop(0)
            l.append(gene)
            if gid in GENOM:
                for item in l:
                    g.write("%s\t" % item)
                g.write("\n")
    f.close()
g.close()

