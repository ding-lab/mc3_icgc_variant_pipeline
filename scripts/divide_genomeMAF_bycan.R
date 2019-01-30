library(data.table)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


#Testing 
#maf = fread("processed_data/genome_broadbed.gaf.maf",header=F,sep="\t",na.strings="NA",colClasses=list(character=c(17,33,35)))
#smap = fread("/diskmnt/Projects/ICGC_MC3/ID_Mapping/MC3_PCAWG_Matched.ids.v3.txt")

maf = fread(args[1],header=F,sep="\t",na.strings="NA",colClasses=list(character=c(17,33,35)))
smap = fread(args[2])



#Reformat the MAF
fmaf = data.frame(
    "Hugo_Symbol"=maf$V43, 
    "Entrez_Gene_Id"=".",	
    "Center"="ICGC","NCBI_Build"="GRCh37",
    "Chromosome"=maf$V1,
    "Start_Position"=maf$V2,
    "End_Position"=maf$V3,
    "Strand"=maf$V4,
    "Variant_Classification"=maf$V5,
    "Variant_Type"=maf$V6,
    "Reference_Allele"=maf$V7,
    "Tumor_Seq_Allele1"=maf$V8,
    "Tumor_Seq_Allele2"=maf$V9,
    "dbSNP_RS"=maf$V10,
    "dbSNP_Val_Status"="NA",
    "Tumor_Sample_Barcode"=maf$V12,
    "Matched_Norm_Sample_Barcode"=maf$V13,
    "Tumor_Sample_UUID"=maf$V42,
    "Cancer"=maf$V41,
    "HGVSc"=maf$V14,
    "Context"=maf$V15,
    "Callers"=maf$V19,
    "t_ref_count"=maf$V37,
    "t_alt_count"=maf$V38
)

fmafi <- merge(fmaf,smap,by.x="Tumor_Sample_Barcode",by.y="tcga_pcawg_aliquot_id")

fmafi$Tumor_Sample_Barcode = fmafi$mc3_exome_barcode

fmafi <- fmafi[c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1,17,18,19,20,21,22,23,24,25,26)]

cans <- unique(fmafi$Cancer)
wigs <- unique(data.frame("id"=fmafi$Tumor_Sample_Barcode,"path"=paste("/diskmnt/Projects/ICGC_MC3/Coverage_Reduction_Genome_Wig/",fmafi$Tumor_Sample_Barcode,".wig",sep="")))


for(i in cans){
    dat <- fmafi[which(fmafi$Cancer == i),]
    oname = paste("processed_data/SMG/input/",i,".genome.reduced.maf",sep="")
    write.table(dat,oname,quote=F,sep="\t",row.names=F)
    samps <- unique(dat$Tumor_Sample_Barcode)
    wig_list <- wigs[which(wigs$id %in% samps),]
    wname = paste("processed_data/SMG/input/",i,".wig_list",sep="")
    write.table(wig_list,wname,quote=F,sep="\t",row.names=F,col.names=F)
}


