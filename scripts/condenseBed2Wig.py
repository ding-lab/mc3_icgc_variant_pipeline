import sys


if __name__ == "__main__":


    #HANDLE THE FIRST CASE
    c,start,end,k,cov = sys.stdin.readline()[:-1].split("\t")
    curr = [c,start,end,k,cov]

    #POSITIONS 
    #0 = CHROMOSOME 
    #1 = START 
    #2 = END
    #3 = ID
    #4 = COVERAGE 


    #HANDLE MIDDLE CASES
    for line in sys.stdin:
        c,start,end,k,cov = line[:-1].split("\t")
        #THIS CATCHES COVERED CHANGES AND SEQUENCE GAPS
        if curr[4] == cov and curr[2] == start:
            curr[2] = end
            curr[3] = k  
        else:
            sys.stdout.write("\t".join(curr))
            sys.stdout.write("\n")
            curr = [c,start,end,k,cov]

    #PRINT THE LAST LINE
    sys.stdout.write("\t".join(curr))
    sys.stdout.write("\n")
