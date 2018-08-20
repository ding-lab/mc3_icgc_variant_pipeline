
library(data.table)
library(UpSetR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

dat = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA")
#This needs a deduplication step 
dat$GENOME_VAF = as.numeric(dat$GENOME_VAF)

#dat2 <- dat[which(dat$FILTER == "PASS" | is.na(dat$FILTER)),] #Depreciated
#mc3Samples <- unique(dat2[which(dat2$FILTER == "PASS"),]$Barcode_MC3)

mc3Samples <- unique(dat$Barcode_MC3)


pdf(file=args[2],height=11,width=22,useDingbat=F)
upset(dat, sets = c("CONCORDANCE","MUTECT","EXOME_MUSE","VARSCAN","RADIA","SOMATICSNIPER","INDELOCATOR","PINDEL","GENOME_MUSE","DKFZ","BROAD","SANGER","SMUFIN"), sets.bar.color = "#BE312D",order.by = "freq", empty.intersections = "on",point.size = 5, line.size = 2,keep.order = T)
dev.off()
             

#Run this one in TMUX because it will likely take forever. 
pdf(file=args[3],height=11,width=72,useDingbat=F)
upset(dat, sets = c("INDELOCATOR","EXOME_MUSE","MUTECT","PINDEL","RADIA","SOMATICSNIPER","VARSCAN","BROAD","DKFZ","GENOME_MUSE","SANGER","SMUFIN","CONCORDANCE"), sets.bar.color = "#BE312D",order.by = "freq", empty.intersections = "on",point.size = 5, line.size = 2,nintersects=NA)
dev.off()

