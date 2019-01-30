

library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(tidyr)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))


data$match = ifelse(!is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$PCAWG_only = ifelse(is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$MC3_only = ifelse(is.na(data$"Chromosome:1") & !is.na(data$Chromosome),1,0)

uu_pcawg <- fread(args[2])
#uu_mc3 <- data[which(data$MC3_only == 1 & (data$FILTER == "PASS" | data$FILTER == "wga")),]
#uu_match <- data[which(data$match == 1 & (data$FILTER == "PASS" | data$FILTER == "wga" | is.na(data$FILTER))),]

uu_mc3 <- data[which(data$MC3_only == 1),]
uu_match <- data[which(data$match == 1),]


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


SNP <- together[which(together$Start == together$End),]

CADD <- fread(args[3],header=T,sep="\t")

all_sc <- merge(SNP,CADD, by.x=c("Chr","Start","Ref","Alt"), by.y=c("CHR", "POS", "REF", "ALT"))

cols <- colnames(all_sc)
facs <- NULL
for(i in cols){
    myc <- class(all_sc[[i]])
    if(myc == "numeric"){
        facs = rbind(facs,i)
    }
}
facs <- as.character(facs[,1])

for(f in facs){
    all_sc$my_bin <- cut(all_sc[[f]],breaks=20,lables=T)
    look4 <- data.frame(all_sc %>% group_by(Uniq,my_bin) %>% summarise (n = n()) %>% mutate(freq = n / sum(n)))

    ft <- ggplot(look4,aes(x=factor(my_bin),y=freq,color=Uniq))
    ft <- ft+geom_point()
    ft <- ft+theme_classic()
    ft <- ft+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position=c(.8,.8))
    ft <- ft+scale_color_manual(values=c("purple","blue","red"))
    ft <- ft+ylab("Uniq Fraction")+xlab(paste("Equally cut ",f," content bins",sep=""))
    ft

    yout <- paste(args[4],f,"_PCAWG_unique.pdf",sep="")
    pdf(yout,height=6,width=6,useDingbats=F)
    print(ft)
    dev.off()
}
