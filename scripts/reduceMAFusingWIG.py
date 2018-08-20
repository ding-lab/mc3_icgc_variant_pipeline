

def makeWIGdict_condensed(wig):
    OUTdict = {}
    ss = 0 
    se = 0 
    with open(wig,"r") as f:
        for line in f:
            chromo,start,end,_,cov = line[:-1].split("\t")
            ss = int(start)+1
            se = int(end)+1
            cov = float(cov)
            if cov == 1:
                if chromo in OUTdict:
                    OUTdict[chromo] += [range(ss,se)]
                else: #first instance of chromo
                    OUTdict[chromo] = [range(ss,se)]
    return OUTdict



def makeWIGdict_oneBase(wig): #Could speed this up by putting in a tree structure... 
    OUTdict = {}
    stretch = False
    ss = 0
    se = 0
    
    with open(wig,"r") as f:
        for line in f:
            chromo,start,end,_,cov = line[:-1].split("\t")
            start = int(start)+1
            end = int(end)+1
            cov = float(cov)
            if cov == 1: # if covered see if in a stretch
                if stretch: # if in coverered region extend
                    se = end
                else: #else set the start
                    ss = start
                stretch = True
            else:
                if stretch: # if prev was covered store range
                    if chromo in OUTdict:
                        OUTdict[chromo] += [range(ss,se)]
                    else: #first instance of chromo
                        OUTdict[chromo] = [range(ss,se)]
                stretch = False
#NOTE no base case. This function isn't called.
 
    return OUTdict


def reduceMAFwithWIG(maf,wig,ome):
    WIG = {}
#    if ome == "exome":
#        WIG = makeWIGdict_oneBase(wig)
    
    if ome == "genome" or ome == "exome":
        WIG = makeWIGdict_condensed(wig)

    with open(maf,"r") as f:
        for line in f:
            chromo,start,end,*_ = line[:-1].split("\t")
            start = int(start)
            end = int(end)
            if chromo in WIG:
                if start == end:
                    if any([start in i for i in WIG[chromo]]):
                        print(line[:-1])
                else: #This is a longer indel and needs to be 50% covered, COULD Change
                    full = end-start+1 #again make this zero based
                    bases_covd = 0
                    for j in range(start,end+1): #0 based
                        bases_covd += sum([j in i for i in WIG[chromo]])
                    if bases_covd >= .5*full: #recipricol overlap 50%
                        print(line[:-1])
                        
                        

if __name__ == "__main__":
    import sys
    maf = sys.argv[1]
    wig = sys.argv[2]
    ome = sys.argv[3] # This is a string of "exome" or "genome". 
    reduceMAFwithWIG(maf,wig,ome)

