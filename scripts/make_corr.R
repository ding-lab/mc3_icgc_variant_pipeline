library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

#Functions 
'%!in%' <- function(x,y)!('%in%'(x,y))


#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}
 
#args = c('processed_data/data.4.cancerOLAP.txt','id_mapping/PCA.sample.cancer.txt','data/MUTperMB.immunepaper.txt','data/PCA.characterization.removed.samples.txt','figures/CancerXConcordance.pdf')


dat = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA")
dat$GENOME_VAF = as.numeric(dat$GENOME_VAF)

dat3 <- dat

code <- fread(args[2],header=F)

dat3$char12 = substr(dat3$Barcode_MC3,1,12)


dat3$char12 = substr(dat3$Barcode_MC3,1,12)


dat4 <- merge(dat3,code,by.x="char12",by.y="V1")

dat4$MATCH <- as.numeric(dat4$CONCORDANCE)
dat4$PCAWG_only <- ifelse(is.na(dat4$FILTER) & dat4$CONCORDANCE == 0, 1, 0)
dat4$MC3_only <- ifelse(!is.na(dat4$FILTER) & dat4$CONCORDANCE == 0, 1, 0)

newdat = data.frame(dat4 %>% group_by(Barcode_MC3) %>% summarise(id_match = sum(CONCORDANCE), id_mc3=sum(MC3_only) , id_pcawg = sum(PCAWG_only)))

###NOW BRING UP GGPLOT CODE 
percent = c("0%","25%","50%","75%","100%")

toPlot <- newdat[which(newdat$id_match+newdat$id_mc3 != 0),]

toPlot$perc_Matched_inall_MC3_vars = toPlot$id_match/(toPlot$id_match+toPlot$id_mc3)
toPlot$perc_Matched_inall_PCAWG_vars = toPlot$id_match/(toPlot$id_match+toPlot$id_pcawg)
toPlot$char12 = substr(toPlot$Barcode_MC3,1,12)
toPlot <- merge(toPlot,code,by.x="char12",by.y="V1")
toPlot2 <- data.frame(toPlot %>% group_by(V2) %>% mutate(n = n())  %>% mutate(label = paste0(V2,"\nN = ",n)))


melt2plot <- melt(toPlot2)
kvar <- c("perc_Matched_inall_MC3_vars","perc_Matched_inall_PCAWG_vars")
myVars2plot <- melt2plot[which(melt2plot$variable %in% kvar),]


p <- ggplot(myVars2plot,aes(label,value,fill=variable))
p <- p+geom_boxplot()
p <- p+theme_classic()
p <- p+theme(legend.position="bottom")
p <- p+scale_fill_brewer("Set1")
p <- p+scale_y_continuous(labels = percent)
p

pdf(args[3],height=4,width=11,useDingbats=F)
print(p)
dev.off()

#This is the correlations part 
#This is the original pass. I need to do another with removal of hypermutators. 

datim <- fread(args[3],sep="\t",header=T)
names(datim) <- gsub("-","_",names(datim))
names(datim) <- gsub(" ","_",names(datim))

rems <- fread(args[4],header=F)
colnames(rems)= c("TCGA_Barcode")
rems$char12 <- substr(rems$TCGA_Barcode,1,12)




imv2p <- merge(myVars2plot,datim,by.x="char12",by.y="Patient_ID")
imv2p_cleaned <- imv2p[which(imv2p$char12 %!in% rems$char12),]


#Raw data 
pmi <- imv2p[which(imv2p$variable=="perc_Matched_inall_PCAWG_vars"),]
pmm <- imv2p[which(imv2p$variable=="perc_Matched_inall_MC3_vars"),]

bycani <- as.data.frame(pmi %>% group_by(V2) %>% summarize(meanPmatchI = mean(value), meanMutI = mean(Non_silent_per_Mb)))
icorr <- cor.test(bycani$meanPmatchI,bycani$meanMutI)
iwilc <- wilcox.test(bycani$meanPmatchI,bycani$meanMutI)
ggplot(bycani, aes(x = meanMutI, y = meanPmatchI)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')


bycanm <- as.data.frame(pmm %>% group_by(V2) %>% summarize(meanPmatchM = mean(value), meanMutM = mean(Non_silent_per_Mb)))
mcorr <- cor.test(bycanm$meanPmatchM,bycanm$meanMutM)
mwilc <- wilcox.test(bycanm$meanPmatchM,bycanm$meanMutM)
ggplot(bycanm, aes(x = meanMutM, y = meanPmatchM)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')

ggplot(bycanm, aes(x = meanMutM, y = meanPmatchM)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm') + 
  ggtitle("739 x 739 MC3 comparison")+
  geom_text(aes(label=V2)) + xlab("Non-Syn Muts per MB") + ylab("MC3 concordance")



#Cleaned out Hypermutators and purity type samples 
pmi <- imv2p_cleaned[which(imv2p_cleaned$variable=="perc_Matched_inall_PCAWG_vars"),]
pmm <- imv2p_cleaned[which(imv2p_cleaned$variable=="perc_Matched_inall_MC3_vars"),]

bycani <- as.data.frame(pmi %>% group_by(V2) %>% summarize(meanPmatchI = mean(value), meanMutI = mean(Non_silent_per_Mb)))
icorr <- cor.test(bycani$meanPmatchI,bycani$meanMutI)
iwilc <- wilcox.test(bycani$meanPmatchI,bycani$meanMutI)
ggplot(bycani, aes(x = meanMutI, y = meanPmatchI)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')


bycanm <- as.data.frame(pmm %>% group_by(V2) %>% summarize(meanPmatchM = mean(value), meanMutM = mean(Non_silent_per_Mb)))
mcorr <- cor.test(bycanm$meanPmatchM,bycanm$meanMutM)
mwilc <- wilcox.test(bycanm$meanPmatchM,bycanm$meanMutM)
ggplot(bycanm, aes(x = meanMutM, y = meanPmatchM)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')


#Now using the rank order of all the TCGA samples
rawTCGA <- data.frame( datim %>% group_by(Cohort) %>% summarize(nspMB = mean(Non_silent_per_Mb)))


pmic <- merge(bycani,rawTCGA,by.x="V2",by.y="Cohort")  
pmmc <- merge(bycanm,rawTCGA,by.x="V2",by.y="Cohort")


wilcox.test(pmic$meanPmatchI,pmic$nspMB)
cor.test(pmic$meanPmatchI,pmic$nspMB)
ggplot(pmic, aes(x = nspMB, y = meanPmatchI)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')


wilcox.test(pmmc$meanPmatchM,pmmc$nspMB)
cor.test(pmmc$meanPmatchM,pmmc$nspMB)
ggplot(pmmc, aes(x = nspMB, y = meanPmatchM)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')


#Now just to check with some of the hypermutators removed 
rawTCGA2 <- data.frame( datim[which(datim$Patient_ID %!in% rems$char12),] %>% group_by(Cohort) %>% summarize(nspMB = mean(Non_silent_per_Mb)))

pmic <- merge(bycani,rawTCGA2,by.x="V2",by.y="Cohort")
pmmc <- merge(bycanm,rawTCGA2,by.x="V2",by.y="Cohort")


pcawgWR = wilcox.test(pmic$meanPmatchI,pmic$nspMB)
pcawgPC = cor.test(pmic$meanPmatchI,pmic$nspMB)
ppcawg <- ggplot(pmic, aes(x = nspMB, y = meanPmatchI)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm')+
  ggtitle("625 x 8852 PCAWG")+
  geom_text(aes(label=V2), hjust = 'left') + xlab("Non-Syn Muts per MB") + ylab("PCAWG concordance") + 
  theme_bw() +
  annotate("text",x=7.5,y=.7,label=paste("Mann-Whitney p-value",signif(pcawgWR$p.value,3),sep="=")) +
  annotate("text",x=7.5,y=.675,label=paste("Pearson p-value",signif(pcawgPC$p.value,3),sep="=")) +  
  annotate("text",x=7.5,y=.65,label=paste("Pearson R^2",signif(pcawgPC$estimate^2,3),sep="="))

pdf("figures/Correlation_PCAWG_MpMB.pdf",height=5.5, width=7)
print(ppcawg)
dev.off()

mc3WR = wilcox.test(pmmc$meanPmatchM,pmmc$nspMB)
mc3PC = cor.test(pmmc$meanPmatchM,pmmc$nspMB)
pmc3 <- ggplot(pmmc, aes(x = nspMB, y = meanPmatchM)) +
  geom_point(shape = 19, size = 2, ) +
  geom_smooth(colour = "red", fill = "lightgreen", method = 'lm') + 
  ggtitle("625 x 8852 MC3")+
  geom_text(aes(label=V2), hjust = 'left') + xlab("Non-Syn Muts per MB") + ylab("MC3 concordance") + 
  theme_bw() +
  annotate("text",x=7.5,y=.7,label=paste("Mann-Whitney p-value",signif(mc3WR$p.value,3),sep="=")) +
  annotate("text",x=7.5,y=.675,label=paste("Pearson p-value",signif(mc3PC$p.value,3),sep="=")) + 
  annotate("text",x=7.5,y=.65,label=paste("Pearson R^2",signif(mc3PC$estimate^2,3),sep="="))

pdf("figures/Correlation_MC3_MpMB.pdf",height=5.5, width=7)
print(pmc3)
dev.off()















