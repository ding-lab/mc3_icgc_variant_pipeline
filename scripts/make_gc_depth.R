library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


#args <- c("removed.by.coverage.genome.maf","CADD.annotations.vars.rem_by_covg.txt","../output/full_cleaned.tsv","unique_icgc.notin.mc3_controlled.txt","CADD.annotations.fullvars.txt")
 
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("insufficient number of input files to R script (input file).n", call.=FALSE)
}


uncovg <- fread(args[1], header=FALSE, na.strings="NA",colClasses=list(character=c(17,33,35)))

uncadd <- fread(args[2], header=T, na.string="NA", colClasses=list(character=c(9)))

full <- fread(input=args[3], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

full$match = ifelse(!is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$PCAWG_only = ifelse(is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$MC3_only = ifelse(is.na(full$"Chromosome:1") & !is.na(full$Chromosome),1,0)


uu_pcawg <- fread(args[4],header=F,na.string="NA", colClasses=list(character=c(135)))
uu_mc3 <- full[which(full$MC3_only == 1),]
uu_match <- full[which(full$match == 1),]

uu_mc3$match <- NULL
uu_mc3$PCAWG_only <- NULL
uu_mc3$MC3_only <- NULL

uu_match$match <- NULL
uu_match$PCAWG_only <- NULL
uu_match$MC3_only <- NULL

uu_pcawg$UNIQUE = "PCAWG"
uu_mc3$UNIQUE = "MC3"
uu_match$UNIQUE = "MATCH"
colnames(uu_pcawg) <- colnames(uu_mc3)

all <- rbind(uu_mc3,uu_pcawg,uu_match)

together = data.frame("Chr"= ifelse(!is.na(all$Chromosome),all$Chromosome,all$"Chromosome:1"))
together$Start <- ifelse(!is.na(all$Start_Position),all$Start_Position,all$"Start_position:1")
together$End <- ifelse(!is.na(all$End_Position),all$End_Position,all$"End_position:1")
together$Ref <- ifelse(!is.na(all$Reference_Allele),all$Reference_Allele,all$"Reference_Allele:1")
together$Alt <- ifelse(!is.na(all$Tumor_Seq_Allele2),all$Tumor_Seq_Allele2,all$"Tumor_Seq_Allele2:1")
together$ID <- all$mc3_exome_barcode
together$Uniq <- all$UNIQUE
together$Depth <- ifelse(!is.na(all$t_depth),all$t_depth,all$"t_alt_count:1"+ all$"t_ref_count:1")
#default to MC3 depth (likely deeper)

SNP <- together[which(together$Start == together$End),]

CADD <- fread(args[5], header=T, sep="\t", na.string="NA", colClasses=list(character=c(9)))

all_sc <- merge(SNP,CADD, by.x=c("Chr","Start","Ref","Alt"), by.y=c("CHR", "POS", "REF", "ALT"))

#This is going to add the mc3 id instead of the PCAWG sample id 
pcawg_mc3 <- unique(data.frame("PCAWG"=full$"Tumor_Sample_Barcode:1","MC3"=full$mc3_exome_barcode))

uncovg_id <- merge(uncovg,pcawg_mc3,by.x="V12",by.y="PCAWG")

un <- data.frame("Chr"=uncovg_id$V1,"Start"=uncovg_id$V2,"Ref"=uncovg_id$V7,"Alt"=uncovg_id$V9,"End"=uncovg_id$V3,"ID"=uncovg_id$MC3,"Uniq"="Removed","Depth"=uncovg_id$V37+uncovg_id$V38) 


all_un <- merge(un,uncadd, by.x=c("Chr","Start","Ref","Alt"), by.y=c("CHR", "POS", "REF", "ALT"))

all <- rbind(all_sc,all_un)
all$my_bin <- cut(all$GC,breaks=20,lables=T)
look4 <- data.frame(all %>% group_by(Uniq,my_bin) %>% summarise(pheno=median(Depth,na.rm=T)))

#Binned data
    ft <- ggplot(look4,aes(x=factor(my_bin),y=pheno,color=Uniq))
    ft <- ft+geom_point()
    ft <- ft+theme_classic()
    ft <- ft+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position=c(.8,.8))
    ft <- ft+scale_color_manual(values=c("purple","blue","red","orange"))
    ft <- ft+ylab("Sequence Depth")+xlab(paste("Equally cut GC content bins",sep=""))
    ft <- ft+ylim(0,150)
    ft

pdf(args[6],height=3,width=4,useDingbats=F)
print(ft)
dev.off()


#Density plot
    p <- ggplot(all,aes(x=GC, fill=Uniq))
    p <- p + geom_density(alpha=0.45)
    p <- p + scale_fill_manual(values=c("purple","blue","red","orange"))
    p <- p + ylab("Density")+xlab("GC")
    p <- p + theme_classic()
    p

pdf(args[7],height=2,width=4,useDingbats=F)
print(p)
dev.off()



#So here are some extra analysis, I'm going to try and capture all of the "removed" variants with Higher GC -content that were "Uncapturable from MC3"

#Here is the code 
c299 <- fread("/diskmnt/Projects/ICGC_MC3/Data/cancer.gene.299.txt",header=F)
un_can <- all_un[which(all_un$GeneName %in% c299$V1),]
table(un_can$ConsDetail)

miss <- un_can[which(un_can$ConsDetail == "missense"),]
extreme_miss <- miss[which(miss$GC < .3 | miss$GC > .7),]

utr3 <- un_can[which(un_can$ConsDetail == "3_prime_UTR"),]
sort(table(utr3$GeneName))


#Now I'm going to attach cancer types because I think this PTPRQ genes is going to be a hit. 
#Phosphatidylinositol phosphatase required for auditory function. May act by regulating the level of phosphatidylinositol 4,5-bisphosphate (PIP2) level in the basal region of hair bundles. Can dephosphorylate a broad range of phosphatidylinositol phosphates, including phosphatidylinositol 3,4,5-trisphosphate and most phosphatidylinositol monophosphates and diphosphates. Phosphate can be hydrolyzed from the D3 and D5 positions in the inositol ring. Has low tyrosine-protein phosphatase activity; however, the relevance of such activity in vivo is unclear. Plays an important role in adipogenesis of mesenchymal stem cells (MSCs). Regulates the phosphorylation state of AKT1 by suppressing the phosphatidylinositol 3,4,5-trisphosphate (PIP3) level in MSCs and preadipocyte cells.

#THAT IS REALLY FUN And I think that it will be abundant in BRCA

cantypes <- fread("/diskmnt/Projects/ICGC_MC3/Gen_Figures/OlapByCancer_Figure/PCA.sample.cancer.txt",header = F)
colnames(cantypes) <- c("char12","CODE")

all_un$char12 <- substr(all_un$ID,1,12)
all_unc <- merge(all_un,cantypes,by="char12")

#Missense exploration
aumiss <- all_unc[which(all_unc$ConsDetail == "missense"),]
yo2 <- aumiss %>% group_by(Chr,Start,CODE) %>% tally()
yo2[order(yo2$n,decreasing=T),]
yo2 <- aumiss %>% group_by(CODE,GeneName) %>% tally()
yo2[order(yo2$n,decreasing=T),]

wdr87skcm <- unique(aumiss[which(aumiss$GeneName == "WDR87" & aumiss$CODE == "SKCM"),]$ID)
wdfy4skcm <- unique(aumiss[which(aumiss$GeneName == "WDFY4" & aumiss$CODE == "SKCM"),]$ID)

#3'UTR Exploration
au3utr <- all_unc[which(all_unc$ConsDetail == "3_prime_UTR"),]
yo2 <- au3utr %>% group_by(Chr,Start,CODE) %>% tally()
yo2[order(yo2$n,decreasing=T),]
yo2 <- au3utr %>% group_by(CODE,GeneName) %>% tally()
yo2[order(yo2$n,decreasing=T),]



#Check these against the all of the controlled MC3 DATA and see if they were removed too. 



