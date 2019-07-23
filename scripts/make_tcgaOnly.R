#Matthew Bailey
#Started January 30, 2019
#Captures variations in SNP/DNP and indel calls in final overlap. 

library(data.table)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(gplots)

#DEVELOPMENT ONLY######
#args <- c("output/full_cleaned.tsv","/diskmnt/Projects/ICGC_MC3/Data/cancer.gene.299.txt")


#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file) of Snakefile rule: snp_tnp_indel_figure ", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

data$match = ifelse(!is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$PCAWG_only = ifelse(is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$MC3_only = ifelse(is.na(data$"Chromosome:1") & !is.na(data$Chromosome),1,0)

can <- fread(args[2],header=F)


to <- data[which(data$MC3_only == 1),]
tof <- to[which(to$FILTER == "PASS"),]
fc <- tof[which(tof$Hugo_Symbol %in% can$V1),]

summary(na.omit(data$t_alt_count/data$t_depth))
summary(na.omit(tof$t_alt_count/tof$t_depth))
summary(na.omit(fc$t_alt_count/fc$t_depth))


badsamp <- c("TCGA-CA-6717-01A-11D-1835-10","TCGA-BR-6452-01A-12D-1800-08")
sum(sort(table(fc[which(!fc$mc3_exome_barcode %in% badsamp ),]$Hugo_Symbol)))
sum(sort(table(fc[which(fc$mc3_exome_barcode %in% badsamp ),]$Hugo_Symbol)))


table(fc$Hugo_Symbol)

pik <- fc[which(fc$Hugo_Symbol == "PIK3CA"),]
kmt2c <- fc[which(fc$Hugo_Symbol == "KMT2C"),]

summary(na.omit(pik$t_alt_count/pik$t_depth))
summary(na.omit(kmt2c$t_alt_count/kmt2c$t_depth))

summary(fc$t_depth)
summary(fc$t_alt)


good <- fc[which(!fc$mc3_exome_barcode %in% badsamp ),]
summary(na.omit(good$t_alt_count/good`$t_depth))
