


def file2set(maff):
    Out = set()
    with open(maff,"r") as f:
        for line in f: 
            Out.add(line)
    f.close()

    return Out



def maf_set_complement(mySet,maff):
    with open(maff,"r") as f:
        for line in f: 
            if line not in mySet:
                print(line[:-1]) #no need for a double '\n'
    f.close()
    

if __name__ == "__main__":
    import sys 
    small_maf = file2set(sys.argv[1])
    large_maf = sys.argv[2]
    maf_set_complement(small_maf,large_maf)
