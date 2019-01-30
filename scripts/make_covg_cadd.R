library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)


#args[1] = pcawg maf removed by coverage 
#args[2] = pcawg rem cadd annotated 
 
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


uncovg <- fread(args[1],header=FALSE, na.strings="NA",colClasses=list(character=c(17,33,35)))
#uncovg <- fread("removed.by.coverage.genome.maf",header=FALSE, na.strings="NA",colClasses=list(character=c(17,33,35)))

#uncadd <- fread()
uncadd <- fread("CADD.annotations.vars.rem_by_covg.txt", header=T, na.string="NA", colClasses=list(character=c(9)))



#full <- fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))
full <- fread("../output/full_cleaned.tsv", sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

full$match = ifelse(!is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$PCAWG_only = ifelse(is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$MC3_only = ifelse(is.na(full$"Chromosome:1") & !is.na(full$Chromosome),1,0)


uu_pcawg <- fread("unique_icgc.notin.mc3_controlled.txt")
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

#CADD <- fread(,header=T,sep="\t")
CADD <- fread("CADD.annotations.fullvars.txt",header=T,sep="\t",na.string="NA", colClasses=list(character=c(9)))

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

#Density plot
    p <- ggplot(all,aes(x=GC, fill=Uniq))
    p <- p + geom_density(alpha=0.45)
    p <- p + scale_fill_manual(values=c("purple","blue","red","orange"))
    p <- p + ylab("Density")+xlab("GC")
    p <- p + theme_minimal()
    p




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
    look4 <- data.frame(all_sc %>% group_by(Uniq,my_bin) %>% summarise(pheno=median(Depth,na.rm=T)))

    ft <- ggplot(look4,aes(x=factor(my_bin),y=pheno,color=Uniq))
    ft <- ft+geom_point()
    ft <- ft+theme_classic()
    ft <- ft+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5),legend.position=c(.8,.8))
    ft <- ft+scale_color_manual(values=c("purple","blue","red"))
    ft <- ft+ylab("Sequence Depth")+xlab(paste("Equally cut ",f," content bins",sep=""))
    ft

    yout <- paste(args[4],f,"_depth.pdf",sep="")
   pdf(yout,height=6,width=6,useDingbats=F)
    print(ft)
    dev.off()
}


