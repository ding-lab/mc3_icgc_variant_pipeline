library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

#For testing purposes
#args = c("output/full_cleaned.tsv","processed_data/Full_Clonality.tsv")

full = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))
 
#Two different ways. 
full$match = ifelse(!is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$PCAWG_only = ifelse(is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$MC3_only = ifelse(is.na(full$"Chromosome:1") & !is.na(full$Chromosome),1,0)


tmp = full[which(full$MC3_only == 1),]
UMC3 = tmp[,c("Chromosome","Start_Position","End_Position","mc3_exome_barcode","Variant_Classification")]
UMC3$uid <- paste(UMC3$Chromosome,UMC3$Start_Position,UMC3$End_Position,UMC3$mc3_exome_barcode,sep="_")
dim(UMC3)

tmp2 = full[which(full$PCAWG_only == 1),]
UPCAWG = tmp2[,c("Chromosome:1","Start_position:1","End_position:1","Donor_ID","Variant_Classification:1")]
UPCAWG$uid <- paste(UPCAWG$"Chromosome:1",UPCAWG$"Start_position:1",UPCAWG$"End_position:1",UPCAWG$"Donor_ID",sep="_")
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

#TODO, and I'm talking to eduard about this one: There are 59,232 mutations without TCGA clonality extimations in these data. Most of these come from a handful of samples but still these were reduced, but hopefully I will be able to recover this for a down stream anlayis., Keep in mind these mutations were restricted to SNV called PCAWG variants. 

#### Now I need to make a similar plot for TCGA/MC3
TCGA_dat <- dat[which(!is.na(dat$"Chromosome")),]
TCGA_dat$TCGA_Uniq = ifelse(is.na(TCGA_dat$"Chromosome:1"),"Unique","Matched")
TCGA_dat$uid <- paste(TCGA_dat$"Chromosome",TCGA_dat$"Start_Position",TCGA_dat$"End_Position",TCGA_dat$mc3_exome_barcode,sep="_")

d4 <- TCGA_dat[,c("uid","clonal.ix")]
d4$Clonal <- ifelse(d4$"clonal.ix"=="TRUE",1,0)
mclone <- d4[which(d4$Clonal == 0),]
UMC3$SubClonal <- ifelse(UMC3$uid %in% mclone$uid,1,0)




