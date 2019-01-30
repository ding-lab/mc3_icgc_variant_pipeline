

def wig2fixedwig(wig,fwig):
    W = "" 
    X = 0 
    Y = 0
    Z = "" 
    
    with open(fwig,"w") as o:
        with open(wig,"r") as f:
            header = f.readline()[:-1]
            W,X,Y,V,Z = f.readline()[:-1].split("\t")
            X = int(X)
            Y = int(Y)
            o.write("track type=wiggle_0 name=SomaticCoverag\n")
            o.write("fixedStep chrom="+W+" start="+str(X)+" step=1\n")
            out = str(int(float(Z)))+"\n"
            for i in range(X,Y):
                o.write(out)
            for line in f: 
                out = ""
                A,B,C,I,D = line[:-1].split("\t")
                B = int(B)
                C = int(C)
                out=str(int(float(D)))+"\n"
                if A == W: #Chromosomes Match
                    if B == Y: #Is incremetal 
                        for i in range(B,C):
                            o.write(out)
                    else: #Not incremental
                        o.write("fixedStep chrom="+A+" start="+str(B)+" step=1\n")
                        for i in range(B,C):
                            o.write(out)
                else: #Not same Chromosome 
                    o.write("fixedStep chrom="+A+" start="+str(B)+" step=1\n")
                    for i in range(B,C):
                        o.write(out)
            #This sets the previous 
                W = A
                X = B
                Y = C
                Z = D
        f.close()
    o.close() 



if __name__ == "__main__":
    import sys 
    wig = sys.argv[1] #This is the Wiggle file that was made from the following 2 steps, bedGraphToBigWig and then bigWig2wig  
    fwig = sys.argv[2] #this is the new output in fixed format, and my step will be 1 for these purposes, I can add this functionality later. 
    wig2fixedwig(wig,fwig)
