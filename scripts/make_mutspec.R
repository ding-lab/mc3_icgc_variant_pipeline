library(data.table)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(gplots)

#DEVELOPMENT ONLY######
#args <- c("processed_data/exome.broadbed.gaf.maf","processed_data/genome_broadbed.gaf.maf", "figures/mc3_mutspec.pdf","figures/pcawg_mutspec.pdf","processed_data/mutspecNotes.txt")


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

mc3 <- fread("processed_data/exome.broadbed.gaf.maf")
mc3 <- mc3[which(!grepl("oxog",mc3$V105)),]


pcawg <- fread("processed_data/genome_broadbed.gaf.maf",header=F,sep="\t",na.strings="NA",colClasses=list(character=c(17,33,35)))
smap = fread("/diskmnt/Projects/ICGC_MC3/ID_Mapping/MC3_PCAWG_Matched.ids.v3.txt")
ctypes = fread("/diskmnt/Projects/ICGC_MC3/ID_Mapping/PCA.sample.cancer.txt",header=F)



#Deal with PCAWG
pmaf <- merge(pcawg,smap,by.x="V12",by.y="tcga_pcawg_aliquot_id")
pmaf$char12 = substr(pmaf$mc3_exome_barcode,1,12)
pmafc <- merge(pmaf,ctypes,by.x="char12",by.y="V1")
mutspsamp <- pmafc %>% group_by(char12) %>% tally()
pmafc$tstv <- paste(pmafc$V7,">",pmafc$V9,sep="")
nonhypers <- mutspsamp[which(mutspsamp$n < 1500),]$char12
pmafc <- pmafc[which(pmafc$char12 %in% nonhypers),]


#DEAL with MC3
mc3$char12 = substr(mc3$V12,1,12)
mc3c <- merge(mc3,ctypes,by.x="char12",by.y="V1")
mutspsamp <- mc3c %>% group_by(V12) %>% tally()
mc3c$tstv <- paste(mc3c$V7,">",mc3c$V9,sep="")
mc3c <- mc3c[which(mc3c$char12 %in% nonhypers),]


#The the tally that I need. 
good = c("A>C","A>G","A>T","C>A","C>G","C>T")
mytstv <- pmafc[which(pmafc$tstv %in% good),]
pcawg_tstv <- data.frame(mytstv %>% group_by(V2.y,tstv) %>% tally())

mytstv <- mc3c[which(mc3c$tstv %in% good),]
mc3_tstv <- data.frame(mytstv %>% group_by(V2.y,tstv) %>% tally())

totest <- merge(mc3_tstv,pcawg_tstv,by=c("V2.y","tstv"))

cancers <- unique(totest$V2.y)

#Cancer type level. 
OUT <- NULL
TOPLOT <- NULL 
for(i in cancers){
    CAN <- totest[which(totest$V2.y == i),]
    CAN$V2.y <- NULL
    row.names(CAN) = CAN$tstv
    CAN$tstv <- NULL
    tCAN <- t(CAN)
    cnt <- rowSums(tCAN)
    tCAN[2,] = tCAN[2,]*(cnt[1]/cnt[2]) #this scales the counts to make sure we differences in shape are not due to sample count differences. 
    out <- chisq.test(tCAN)
#    contrib <- 100*out$residuals^2/out$statistic
    contrib <- 100*out$residuals^2
    pf <- data.frame("Cancer"=i,t(contrib[1,]))
    TOPLOT <- rbind(TOPLOT,pf)
    df <- data.frame("Cancer"=i,"stat"=out$statistic,"pval"=out$p.value)
    OUT <- rbind(OUT,df)
}

row.names(TOPLOT) <- TOPLOT$Cancer 
TOPLOT$Cancer <- NULL 
pdf("figures/balloon.mutSpec.pdf",height=6,width=6,useDingbats=F)
balloonplot(t(as.table(as.matrix(TOPLOT))), main ="ChiSq residuals^2 contribution", xlab ="", ylab="",label = FALSE, show.margins = FALSE)
dev.off()


#Look into a couple cancer typs to get numbers 
i = "KICH"



##############################################################################################################3
#Now I want to do this with the unique calls 
#Dont forget to remove hypermutators 

full <- fread("/diskmnt/Projects/ICGC_MC3/mc3_icgc_variant_pipeline/output/full_cleaned.tsv",sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))
ctypes = fread("/diskmnt/Projects/ICGC_MC3/ID_Mapping/PCA.sample.cancer.txt",header=F)

full$match = ifelse(!is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$PCAWG_only = ifelse(is.na(full$Chromosome) & !is.na(full$"Chromosome:1"),1,0)
full$MC3_only = ifelse(is.na(full$"Chromosome:1") & !is.na(full$Chromosome),1,0)

fullm <- full[which(full$MC3_only == 1 | full$match == 1),]
fullp <- full[which(full$PCAWG_only == 1 | full$match == 1),]

#DEAL with PCAWG 
fullp$char12 = substr(fullp$mc3_exome_barcode,1,12)
pull <- merge(fullp, ctypes, by.x="char12",by.y="V1")
mutspsamp <- pull %>% group_by(char12) %>% tally()
nonhypers <- mutspsamp[which(mutspsamp$n < 1500),]$char12 #This repesnets the 95% ile 
pull$tstv <- paste(pull$"Reference_Allele:1",">",pull$"Tumor_Seq_Allele2:1",sep="")
pull <- pull[which(pull$char12 %in% nonhypers),]

#DEAL with MC3
fullm$char12 = substr(fullm$mc3_exome_barcode,1,12)
mull <- merge(fullm, ctypes, by.x="char12",by.y="V1")
mutspsamp <- mull %>% group_by(char12) %>% tally()
mull$tstv <- paste(mull$"Reference_Allele",">",mull$"Tumor_Seq_Allele2",sep="")
mull <- mull[which(pull$char12 %in% nonhypers),]


good = c("A>C","A>G","A>T","C>A","C>G","C>T")
#get counts for these 
pull_tt = pull[which(pull$tstv %in% good),]
pcnt = data.frame(pull_tt %>% group_by(V2,tstv) %>% tally())

mull_tt = mull[which(mull$tstv %in% good),]
mcnt = data.frame(mull_tt %>% group_by(V2,tstv) %>% tally())

#Pull together 
totull <- merge(mcnt,pcnt,by=c("V2","tstv"),all=T)
totull[is.na(totull)] <- 0


#Cancer type level.
cancers <- unique(totull$V2) 
OUTULL <- NULL
TOPLOTULL <- NULL

for(i in cancers){
    CAN <- totull[which(totull$V2 == i),]
    CAN$V2 <- NULL
    row.names(CAN) = CAN$tstv
    CAN$tstv <- NULL
    tCAN <- t(CAN)
    cnt <- rowSums(tCAN)
    tCAN[2,] = tCAN[2,]*(cnt[1]/cnt[2]) #this scales the counts to make sure we differences in shape are not due to sample count differences. 
    out <- chisq.test(tCAN)
#    contrib <- 100*out$residuals^2/out$statistic
    contrib <- 100*out$residuals^2
    pf <- data.frame("Cancer"=i,t(contrib[1,]))
    TOPLOTULL <- rbind(TOPLOTULL,pf)
    df <- data.frame("Cancer"=i,"stat"=out$statistic,"pval"=out$p.value)
    OUTULL <- rbind(OUTULL,df)
}

@SOMETHING IS GOING WRONG HERE! I DON'T' KNOW WHAT.... 
#####################################################################################################################33


#IT Also looks like there may be some problemsome cancertypes from the information above. 
COAD,KICH,LIHC,LUAD,LUSC,OV,READ,STAD 

#PCAWG sample level
mytstv <- pmafc[which(pmafc$tstv %in% good),]
pcawg_tstv_s <- data.frame(mytstv %>% group_by(char12,V2.y,tstv) %>% tally())

#MC3 samples level 
mytstv <- mc3c[which(mc3c$tstv %in% good),]
mc3_tstv_s <- data.frame(mytstv %>% group_by(char12,V2.y,tstv) %>% tally())

totest_s <-  merge(mc3_tstv_s,pcawg_tstv_s,by=c("char12","V2.y","tstv"),all=T)
otest_s[is.na(totest_s)] <- 0
samps <- unique(totest_s$char12)

SAMPS = NULL
for(i in samps){
    S <- totest_s[which(totest_s$char12 == i),]
    row.names(S) = S$tstv
    code = unique(paste(S$char12,"_",S$V2.y,sep=""))
    S$V2.y <- NULL
    S$tstv <- NULL
    S$char12 <- NULL
    tS <- t(S)
    print(code)
    print(tS)
    out <- chisq.test(tS)
    df <- data.frame("ID"=code,"stat"=out$statistic,"pval"=out$p.value)
    SAMPS <- rbind(SAMPS,df)
}

#This is a list of significant differences and the individual level, 
samps = c("TCGA-A6-2683","TCGA-AF-2689","TCGA-BL-A13J","TCGA-AG-3885","TCGA-49-4486","TCGA-44-2659","TCGA-KN-8418")



