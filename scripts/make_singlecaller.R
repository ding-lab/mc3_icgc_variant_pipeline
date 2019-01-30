library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


mc3 <- fread(args[1],header=F,sep="\t")

nonsingle_mc3 <- mc3[which(grepl("\\|", mc3$V111)),]
single_mc3 <- mc3[which(!grepl("\\|", mc3$V111)),]


#This is set of scripts to capture the information needed to plot the composition of the single caller
filt_sc = data.frame(single_mc3 %>% group_by(V109,V111) %>% tally())

pass_tmp <- filt_sc[which(grepl("PASS",filt_sc$V109)),]
pass <- data.frame(pass_tmp %>% group_by(V111) %>% tally())
rest <- filt_sc[which(!grepl("PASS",filt_sc$V109)),]
pass$Filter = "PASS"


npp_tmp <- rest[which(grepl("nonpreferredpair",rest$V109)),]
npp <- data.frame(npp_tmp %>% group_by(V111) %>% tally() )
rest <- rest[which(!grepl("nonpreferredpair",rest$V109)),]
npp$Filter = "nonpreferredpair"


wga_tmp <- rest[which(grepl("wga", rest$V109)),]
wga <- data.frame(wga_tmp %>% group_by(V111) %>% tally())
rest <- rest[which(!grepl("wga", rest$V109)),]
wga$Filter = "wga"


bpon_tmp <- rest[which(grepl("broad_PoN_v2", rest$V109)),]
bpon <- data.frame(bpon_tmp %>% group_by(V111) %>% tally())
rest <- rest[which(!grepl("broad_PoN_v2", rest$V109)),]
bpon$Filter = "broad_PoN_v2"


oxog_tmp <- rest[which(grepl("oxog", rest$V109)),]
oxog <- data.frame(oxog_tmp %>% group_by(V111) %>% tally())
rest <- rest[which(!grepl("oxog",rest$V109)),] 
oxog$Filter = "oxog"


common_tmp <- rest[which(grepl("common_in_exac",rest$V109)),]
common <- data.frame(common_tmp %>% group_by(V111) %>% tally())
rest <- rest[which(!grepl("common_in_exac",rest$V109)),]
common$Filter = "common_in_exac"


other <- data.frame(rest %>% group_by(V111) %>% tally())
other$Filter = "Other"
rest <- NULL

#These are the last too
#pcadont_tmp <- rest[which(grepl("pcadontuse",rest$V109)),]
#pcadont <- data.frame(pcadont_tmp %>% group_by(V111) %>% tally())
#rest <- rest[which(!grepl("pcadontuse",rest$V109)),]
#pcadont$Filter = "pcadontuse"


#strand_tmp <- rest[which(grepl("StrandBias",rest$V109)),]
#strand <- data.frame(strand_tmp %>% group_by(V111) %>% tally())
#rest <- rest[which(!grepl("StrandBias",rest$V109)),]
#strand$Filter = "StrandBias"


toPlot <- rbind(pass,npp,wga,bpon,oxog,common,other)
toPlot$FilterOrd <- factor(toPlot$Filter,levels=c("PASS","nonpreferredpair","wga","broad_PoN_v2","oxog","common_in_exac","Other"))
toPlot$CallerOrd <- factor(toPlot$V111,levels=c("MUTECT","MUSE","RADIA","VARSCANS","SOMATICSNIPER","INDELOCATOR","PINDEL","VARSCANI"))

s <- ggplot(toPlot,aes(x=CallerOrd,y=nn,fill=FilterOrd))
s <- s+geom_bar(stat="identity",position="dodge",width=.75)
s <- s+theme_classic()
s <- s+geom_vline(xintercept=seq(1.5,8.5,1),color="grey")
s <- s+scale_fill_manual(values=brewer.pal("Paired",n=8))
s <- s+coord_flip()
s <- s+theme(legend.position = c(0.8, 0.6))
s <- s+xlab("")+ylab("sqrt(Count)")
s <- s+scale_y_sqrt(breaks=c(0,32,64,128,256,512,1024,2048,4096))
s


pdf(args[2],height=4,width=6,useDingbats=F)
print(s)
dev.off()


