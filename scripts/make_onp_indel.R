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
#args <- c("output/full_cleaned.tsv")


#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file) of Snakefile rule: snp_tnp_indel_figure ", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))


dat <- data %>% group_by(data$Variant_Type, data$"Variant_Type:1") %>% tally()
colnames(dat) <- c("MC3","ICGC","n")

dat$ID = paste(dat$MC3,dat$ICGC,sep=":")


p <- ggplot(dat,aes(x=reorder(ID, -n),y=log10(n+1)))
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=.5))
p <- p + geom_text(aes(label=n), vjust=0)


pdf(args[2],useDingbats=F,height=4,width=4.5)
print(p)
dev.off()
