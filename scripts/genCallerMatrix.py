import os
import sys


fname = sys.argv[1]


with open(fname,"r") as f:
    out = ["Barcode_MC3","INDELOCATOR","EXOME_MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCAN","BROAD","DKFZ","GENOME_MUSE","SANGER","SMUFIN","CONCORDANCE","EXOME_VAF","GENOME_VAF","FILTER"]
    print("\t".join(out))
    f.readline()
    for line in f:
        l = line.strip().split("\t")
        #ids
        exome = l[110]
        genome = l[153]
        ecall = l[106]
        gcall = l[129] 
        filt = l[104]
        evaf = 0
        gvaf = 0
        myid = l[154]
        if exome != "NA":
            evaf = float(l[37])/float(l[35])
        if genome != "NA":
            gvaf = l[138]
            if "|" in gvaf:
                a = gvaf.strip().split("|")
                s = len(a)
                b = sum([float(i) for i in a])
                gvaf = b/s
        #Callers
        indelocator = 0
        emuse = 0
        mutect = 0
        pindel = 0
        radia = 0
        somaticsniper = 0
        varscan = 0        
        broad = 0
        dkfz = 0
        gmuse = 0
        sanger = 0
        smufin = 0
        concordance = 0

        if "INDELOCATOR" in ecall:
            indelocator = 1

        if "MUSE" in ecall:
            emuse = 1

        if "MUTECT" in ecall:
            mutect = 1

        if "PINDEL" in ecall:
            pindel = 1

        if "RADIA" in ecall:
            radia = 1

        if "SOMATICSNIPER" in ecall:
            somaticsniper = 1

        if "VARSCAN" in ecall:
            varscan = 1

        if "broad" in gcall:
            broad = 1

        if "dkfz" in gcall:
            dkfz = 1

        if "muse" in gcall:
            gmuse = 1

        if "sanger" in gcall:
            sanger = 1

        if "smufin" in gcall:
            smufin = 1
     
        if "NA" in ecall or "NA" in gcall:
            concordance = 0
        else:
            concordance = 1
        
        out = [myid,indelocator,emuse,mutect,pindel,radia,somaticsniper,varscan,broad,dkfz,gmuse,sanger,smufin,concordance,evaf,gvaf,filt]
        print("\t".join(map(str,out)))
f.close()
