

def likertFormatting(olapf):
    myFILT = {} 
    filters = ["common_in_exac","gapfiller","nonpreferredpair","native_wga_mix","wga","StrandBias"] 
    PASS = 0
    NA = 0
    
    with open(olapf,"r") as f:
        for line in f:
            inExome = False
            inGenome = False 
            l = line[:-1].split("\t") 
            filt = l[104] #Column for
            if l[0] != "NA":
                inExome = True
            if l[111] != "NA":
                inGenome = True

            for i in filters:
                if i in filt:
                    if i in myFILT: #myDB of counting for each filter
                        if inGenome: 
                            myFILT[i][1] += 1 #increase genome
                        else:
                            myFILT[i][0] += 1 #increase exome only
                    else:
                        if inGenome:
                            myFILT[i] = [0,1] #set genome first 
                        else:
                            myFILT[i] = [1,0] #set exome first 
                else:
                    if filt == "PASS":
                        PASS += 1
                    if filt == "NA":
                        NA += 1
    f.close()
    
#Now I'm going to print this statement and go from there. 
    
    for j in myFILT:
        for k in range(1,myFILT[j][0]+1):
            out = [j,"1"]
            print("\t".join(out))
        for k in range(1,myFILT[j][1]+1):
            out = [j,"2"]
            print("\t".join(out))

                
if __name__ == "__main__":
   import sys 
   olapf = sys.argv[1]
   likertFormatting(olapf)
