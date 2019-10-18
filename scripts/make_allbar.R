library(data.table)	
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)
library(UpSetR)
library(grid)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

#For testing purposes
#args = c("output/full_cleaned.tsv","processed_data/Full_Clonality.tsv","id_mapping/PCA.sample.cancer.txt","processed_data/CADD.annotations.fullvars.txt")

full = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))
 
#Two different ways. 
full$match = ifelse(!is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$PCAWG_only = ifelse(is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$MC3_only = ifelse(is.na(full$"Chromosome:1") & !is.na(full$Chromosome),1,0)

#GET CANCER TYPES 
cans <- fread(args[3],header=F)
colnames(cans) <- c("char12","Cancer")
badcan <- c("THCA","KICH","PRAD")

#GET CADD ANNOTATIONS
CADD <- fread(args[4], header=T, sep="\t", na.string="NA", colClasses=list(character=c(9)))
CADD$ChromPos <- paste(CADD$CHR,CADD$POS,sep=":") 
CADD2add <- CADD[,c("ChromPos","GC")]
CADDhighlow <- CADD2add[which(CADD2add$GC > .70 | CADD2add$GC < .30),]


tmp = full[which(full$MC3_only == 1),]
tmp$VAF <- tmp$t_alt_count/tmp$t_depth
UMC3 = tmp[,c("Chromosome","Start_Position","End_Position","mc3_exome_barcode","Variant_Classification","VAF","CENTERS")]
UMC3$uid <- paste(UMC3$Chromosome,UMC3$Start_Position,UMC3$End_Position,UMC3$mc3_exome_barcode,sep="_")
UMC3$char12 <- substr(UMC3$mc3_exome_barcode,1,12)
UMC3$ChromPos <- paste(UMC3$Chromosome,UMC3$Start_Position,sep=":")
UMC3 <- merge(UMC3,cans,all.x=T)
dim(UMC3)

tmp2 = full[which(full$PCAWG_only == 1),]
#tmp2$VAF #This is actuall coded as i_VAF
UPCAWG = tmp2[,c("Chromosome:1","Start_position:1","End_position:1","mc3_exome_barcode","Donor_ID","Variant_Classification:1","i_VAF","i_Callers")]
UPCAWG$uid <- paste(UPCAWG$"Chromosome:1",UPCAWG$"Start_position:1",UPCAWG$"End_position:1",UPCAWG$"Donor_ID",sep="_")
UPCAWG$char12 <- substr(UPCAWG$mc3_exome_barcode,1,12)
UPCAWG$ChromPos <- paste(UPCAWG$"Chromosome:1",UPCAWG$"Start_position:1",sep=":")
UPCAWG <- merge(UPCAWG,cans)
dim(UPCAWG)


#Not capture the subclonal PCAWG vars
dat = fread(input=args[2], sep="\t", header=TRUE, na.strings="NA")
ui <- dat[which(is.na(dat$Chromosome)),]
um <- dat[which(is.na(dat$"Chromosome:1")),]
match <- dat[which(!is.na(dat$Chromosome) & !is.na(dat$"Chromosome:1") ),]

match$clusternum <- as.numeric((apply(match[,176:180],1,which.max)))
dat$clusternum <- as.numeric((apply(dat[,176:180],1,which.max)))

#### MAKE TABLE FOR SUBCLONE PCAWG ####
ICGC_dat <- dat[which(!is.na(dat$"Chromosome:1")),]
ICGC_dat$ICGC_Uniq = ifelse(is.na(ICGC_dat$Chromosome),"Unique","Matched")
ICGC_dat$uid <- paste(ICGC_dat$"Chromosome:1",ICGC_dat$"Start_position:1",ICGC_dat$"End_position:1",ICGC_dat$"Donor_ID",sep="_")

d3 <- ICGC_dat[,c("uid","clusternum")]
d3$Clonal <- ifelse(d3$clusternum > 1, 1, 0)
pclone <- d3[which(d3$Clonal == 1),]
#### ASSIGN PCAWG-unique with Clonal ### info 
UPCAWG$SubClonal <- ifelse(UPCAWG$uid %in% pclone$uid,1,0)
UPCAWG$VAF10 <- ifelse(UPCAWG$i_VAF < .10 & UPCAWG$i_VAF >= 0.05,1,0)
UPCAWG$VAF5 <- ifelse(UPCAWG$i_VAF < 0.05,1,0)
p_indels <- c("Frame_Shift_Del","Frame_Shift_Ins","De_novo_Start_InFrame","Start_Codon_Ins","Stop_Codon_Ins","In_Frame_Del","In_Frame_Ins","Stop_Codon_Del","Start_Codon_Del")
p_missense <- c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site")
p_other <- c("5'UTR","RNA","5'Flank","Silent","3'UTR","Intron","IGR","lincRNA","De_novo_Start_OutOfFrame","Start_Codon_SNP")

UPCAWG$Indels <- ifelse(UPCAWG$"Variant_Classification:1" %in% p_indels,1,0)
UPCAWG$MissensePlus <- ifelse(UPCAWG$"Variant_Classification:1" %in% p_missense,1,0)
UPCAWG$Other_VarClass <- ifelse(UPCAWG$"Variant_Classification:1" %in% p_other,1,0)
UPCAWG$MMcomplement <- ifelse(! grepl("muse",UPCAWG$i_Callers) & ! grepl("broad",UPCAWG$i_Callers),1,0)
UPCAWG$THCA_KICH_PRAD <- ifelse(UPCAWG$Cancer %in% badcan,1,0)
UPCAWG$GCcontent <- ifelse(UPCAWG$ChromPos %in% CADDhighlow$ChromPos,1,0)

#FOR PLOTTING!!! 
up2plot <- UPCAWG[,c("uid","SubClonal","VAF10","VAF5","Indels","MissensePlus","Other_VarClass","MMcomplement","THCA_KICH_PRAD","GCcontent")]

pdf("figures/SourcesPCAWG.upsetR.pdf",useDingbats=F,height=6,width=12)
upset(up2plot,order.by = "freq", sets = c("SubClonal","VAF10","VAF5","Indels","MissensePlus","Other_VarClass","MMcomplement","THCA_KICH_PRAD","GCcontent"))
grid.text("Sources of variation for PCAWG unique calls",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()



#TODO, and I'm talking to eduard about this one: There are 59,232 mutations without TCGA clonality extimations in these data. Most of these come from a handful of samples but still these were reduced, but hopefully I will be able to recover this for a down stream anlayis., Keep in mind these mutations were restricted to SNV called PCAWG variants. 

#### Now I need to make a similar plot for TCGA/MC3
TCGA_dat <- dat[which(!is.na(dat$"Chromosome")),]
TCGA_dat$TCGA_Uniq = ifelse(is.na(TCGA_dat$"Chromosome:1"),"Unique","Matched")
TCGA_dat$uid <- paste(TCGA_dat$"Chromosome",TCGA_dat$"Start_Position",TCGA_dat$"End_Position",TCGA_dat$mc3_exome_barcode,sep="_")

d4 <- TCGA_dat[,c("uid","clonal.ix")]
d4$Clonal <- ifelse(d4$"clonal.ix"=="TRUE",1,0)
mclone <- d4[which(d4$Clonal == 0),]
UMC3$SubClonal <- ifelse(UMC3$uid %in% mclone$uid,1,0)
UMC3$VAF10 <- ifelse(UMC3$VAF < 0.1 & UMC3$VAF >= 0.05,1,0)
UMC3$VAF5 <- ifelse(UMC3$VAF < 0.05, 1, 0)

m_indels <-  c("Frame_Shift_Del","In_Frame_Ins","Frame_Shift_Ins","In_Frame_Del")
m_missense <- c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site")
m_other <- c("RNA","3'UTR","5'UTR","5'Flank","Silent","3'Flank","Intron","Translation_Start_Site")

UMC3$Indels <- ifelse(UMC3$Variant_Classification %in% m_indels,1,0)
UMC3$MissensePlus <- ifelse(UMC3$Variant_Classification %in% m_missense,1,0)
UMC3$Other_VarClass <- ifelse(UMC3$Variant_Classification %in% m_other,1,0) 
UMC3$MMcomplement <- ifelse(!grepl("MUSE",UMC3$CENTERS) & !grepl("MUTECT",UMC3$CENTERS),1,0)
UMC3$THCA_KICH_PRAD <- ifelse(UMC3$Cancer %in% badcan,1,0)
UMC3$GCcontent <- ifelse(UMC3$ChromPos %in% CADDhighlow$ChromPos,1,0)

#THIS IS FOR THE PLOTTING 
um2plot <- UMC3[,c("uid","SubClonal","VAF10","VAF5","Indels","MissensePlus","Other_VarClass","MMcomplement","THCA_KICH_PRAD","GCcontent")]

pdf("figures/SourcesMC3.upsetR.pdf",useDingbats=F,height=6,width=12)
upset(um2plot,order.by = "freq", sets = c("SubClonal","VAF10","VAF5","Indels","MissensePlus","Other_VarClass","MMcomplement","THCA_KICH_PRAD","GCcontent"))
grid.text("Sources of variation for MC3 unique calls",x = 0.65, y=0.95, gp=gpar(fontsize=20))
dev.off()

