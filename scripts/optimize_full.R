library(data.table)
require(reshape2)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA")
#This needs a deduplication step 



samples = unique(data$mc3_exome_barcode)
samples = unique(substr(data$mc3_exome_barcode,1,12))

CLEAN <- data[which(!grepl("oxog",data$FILTER)),]
CLEAN$match = ifelse(!is.na(CLEAN$Chromosome) & !is.na(CLEAN$"Chromosome:1"),1,0)
CLEAN$PCAWG_only = ifelse(is.na(CLEAN$Chromosome) & !is.na(CLEAN$"Chromosome:1"),1,0)
CLEAN$MC3_only = ifelse(is.na(CLEAN$"Chromosome:1") & !is.na(CLEAN$Chromosome),1,0)

PASS <- CLEAN[which(!duplicated(CLEAN[,c(1:3,12)]) | !duplicated(CLEAN[,c(112:154,156:158)])),]
#this puts the variants in order
oPASS <- PASS[order(PASS$Chromosome,PASS$Start_Position,PASS$End_Position,PASS$"Chromosome:1",PASS$"Start_position:1",PASS$"End_position:1"),]

#This captures the PASS events that are up duplicated in PCAWG and unique in MC3 DNP events
dpPASS <- oPASS[which(!duplicated(oPASS[,c(1:3,12)]) & (duplicated(oPASS[,c(112:154,156:158)]) | duplicated(oPASS[,c(112:154,156:158)],fromLast=TRUE)) & !is.na(oPASS$"Chromosome:1") ),]

#This captures the MC3 events that are duplicated in MC3 and unique in PCAWG
dmPASS <- oPASS[which( (duplicated(oPASS[,c(1:3,12)]) | duplicated(oPASS[,c(1:3,12)],fromLast=TRUE) ) & !duplicated(oPASS[,c(112:154,156:158)])  & !is.na(oPASS$Chromosome) ),]

#Then we can remove these events, The strategy here is to take a single even within the larger events for comparison 
 dp1event <- dpPASS[cumsum(rle(paste(dpPASS$"Start_position:1",dpPASS$"Donor_ID",sep="_"))$lengths),]

 dm1event <- dmPASS[cumsum(rle(paste(dmPASS$Start_Position,dmPASS$Tumor_Sample_Barcode,sep="_"))$lengths),]


#And now I can take out the dmPASS and dpPASS and insert a single event 
o_dp_PASS <- setdiff(oPASS,dpPASS)
o_dp_dm_PASS <- setdiff(o_dp_PASS,dmPASS)

P1 <- rbind(o_dp_dm_PASS,dp1event)
P2 <- rbind(P1,dm1event)

P2$match <- NULL
P2$PCAWG_only <- NULL
P2$MC3_only <- NULL

write.table(P2,file=args[2],sep="\t",quote=F,row.names=F)

