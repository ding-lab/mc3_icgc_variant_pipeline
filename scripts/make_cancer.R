library(data.table)
library(UpSetR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)


#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

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


