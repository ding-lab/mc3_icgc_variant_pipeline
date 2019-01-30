def get_one_headerinfo(afname):
# This function takes a VCF and amkes an ordered header to keep consistency across the files.
    extra_head = False
    out = []
    with open(afname, "r") as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("## Extra column keys:"):
                    extra_head = True
                    continue
                if extra_head:
                     out.append(line.split(" ",2)[1])
            else:
                break
    f.close()

    pout = ["CHR","POS","REF","ALT"]
    for o in out:
        pout.append(str(o))
    print("\t".join(pout))
    return out 


def reformat_anno(afname,ordExtra):
#This needs an order:
    extra = False
    Vars = {} 
    with open(afname, "r") as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("## Extra column keys:"):
                    extra = True
                if extra:
                    continue
            elif line.startswith("#Uploaded_variation"):
                extra = False
                header = line[:-1].split("\t")
            else:
                Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,E = line[:-1].split("\t")
                chrom,pos,ra = Uploaded_variation.split("_")
                ref,alt = ra.split("/")
                ex = E.split(";")
                eHash = {}
                for val in ex:
                    k,v = val.split("=")
                    eHash[k] = v
                annoline = [chrom,pos,ref,alt]
                for e in ordExtra:
                    if e in eHash:
                        annoline.append(str(eHash[e]))
                    else:
                        annoline.append("NA")
                print("\t".join(annoline))
    f.close()




if __name__ == "__main__":
    import sys 
    anno1 = sys.argv[1]
    Extra = get_one_headerinfo(anno1)
    anno2 = sys.argv[2]
    anno3 = sys.argv[3]
    annos = [anno1,anno2,anno3]
    for a in annos:
        reformat_anno(a,Extra)
