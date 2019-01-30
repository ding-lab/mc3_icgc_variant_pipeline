#Matthew Bailey#
#May 30th, 2018#
#Updated Aug 16, 2018
#Python 3# 
#Capture controlled calls made by mc3 in GDC uniq calls# 


def file2list(cfname,h):
    out = []
    with open(cfname, "r") as f:
        if h:
            cfhead = f.readline()
        for line in f:
            out.append(line[:-1])
    return out


def CheckMAFinMAF(full,controlled,MC3inPCAWG,PCAWGnotinMC3,PCAWGinMC3):
    uniq_maf1 = {} 
    uniq_maf2 = {}

    controlled_reported = [] 

    with open(full,"r") as f:
        ln=0
        for line in f:
            ln += 1
            l = line[:-1].split("\t")
            # First i want to check if this is missing in MC3 or GDC 
            mc3filt = l[104]
#            if "oxog" in mc3filt:
#                next
            #MC3
            m_chr = l[0]
            m_start = l[1]
            m_stop = l[2]
            m_id = l[11]
            
#            check = [m_chr,m_start,m_stop,m_id]
#            print(check)

            #GDC 
            g_chr = l[111] #this does have chr
            g_start = l[112]
            g_stop = l[113]
            g_id = l[154]
#            check = [g_chr,g_start,g_stop,g_id]
#            print(check)
#            if ln == 4:
#                sys.exit()

            # Next I need to pick something apart to see if there is anything different on 
            if m_chr == "NA":
                a = [g_chr,int(g_start),int(g_stop),line,""] #this doesn't have 'chr'
                if g_id in uniq_maf2:
                    uniq_maf2[g_id].append(a)
                else:
                    uniq_maf2[g_id] = [a]
            elif g_chr == "NwA":
                b = [m_chr,int(m_start),int(m_stop),line,""] #this doesn't have 'chr'
                if m_id in uniq_maf1:
                    uniq_maf1[m_id].append(b)
                else:
                    uniq_maf1[m_id] = [b]
#            else:
 #               print(line[:-1]) #This prints the full line if the CALL was seen by MC3, thus capturing MATCH and uniq to MC3
    f.close() 

    with open(controlled, "r") as f:
        for line in f:
            l = line[:-1].split("\t")
            c_chr = l[4]
            c_start = l[5]
            c_stop = l[6]
            c_id = l[15]

            if c_id in uniq_maf2:
                for i in uniq_maf2[c_id]:
                    if i[0] == c_chr:
                        if int(c_start) == i[1]:
                            i[4] = line
                            controlled_reported.append(i)
         #                   print(l)
                        elif int(c_stop) == i[2]:
                            i[4] = line
                            controlled_reported.append(i)
         #                   print(l)
                        elif int(c_stop) >= i[1] and int(c_start) < i[2]:
                            i[4] = line
                            controlled_reported.append(i)
        #                    print(l) #All of these print(l) are to capture the controlled MAF
    f.close()


#This is a write out for opposite MAF 
    with open(MC3inPCAWG,"w") as mc:
        with open(PCAWGnotinMC3,"w") as no:
            with open(PCAWGinMC3,"w") as o:
                for i in uniq_maf2:
                    for j in uniq_maf2[i]:
                        if j in controlled_reported:
                            o.write(j[3])
                            mc.write(j[4])
                        else:
                            no.write(j[3]) #Actually Unique to MAF2
            o.close()
        no.close()
    mc.close()


if __name__ == "__main__":
    import sys
    import gzip
    fullMAFs = sys.argv[1]
    controlled = sys.argv[2]
    MC3inPCAWG=sys.argv[3]
    PCAWGnotinMC3=sys.argv[4]
    PCAWGinMC3=sys.argv[5]

    CheckMAFinMAF(fullMAFs,controlled,MC3inPCAWG,PCAWGnotinMC3,PCAWGinMC3)


