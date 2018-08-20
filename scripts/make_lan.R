library(data.table)
require(reshape2)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(scales)
library(dplyr)
library(stringr)
library(RColorBrewer)

#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

#Two different ways. 
data$match = ifelse(!is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$PCAWG_only = ifelse(is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$MC3_only = ifelse(is.na(data$"Chromosome:1") & !is.na(data$Chromosome),1,0)

samples = unique(data$mc3_exome_barcode)
samples = unique(substr(data$mc3_exome_barcode,1,12))

newdat = data.frame(data %>% group_by(mc3_exome_barcode) %>% summarise(id_match = sum(match), id_mc3=sum(MC3_only) , id_pcawg = sum(PCAWG_only)))
#plot(y=(newdat$id_match/(newdat$id_match+newdat$id_mc3)),x=(newdat$id_match/(newdat$id_match+newdat$id_pcawg)))


#From the filter section pick the strategy that works best for this particular group and go from there. This is now being taken care of in the output file 
#PASS_P <- data[which(grepl("PASS",data$FILTER) | is.na(data$FILTER)),]
#PASS_O <- data[which(!grepl("oxog",data$FILTER)),]
#PASS_E = PASS_O

PASS = data

newdat = data.frame(PASS %>% group_by(mc3_exome_barcode) %>% summarise(id_match = sum(match), id_mc3=sum(MC3_only) , id_pcawg = sum(PCAWG_only)))
#plot(y=(newdat$id_match/(newdat$id_match+newdat$id_mc3)),x=(newdat$id_match/(newdat$id_match+newdat$id_pcawg)))

#plot(y=(newdat$id_match/(newdat$id_match+newdat$id_mc3)),x=(newdat$id_match/(newdat$id_match+newdat$id_pcawg)))
###NOW BRING UP GGPLOT CODE 
percent = c("0%","25%","50%","75%","100%")

toPlot <- newdat[which(newdat$id_match+newdat$id_mc3 != 0),]

toPlot$perc_PCAWG_in_MC3 = toPlot$id_match/(toPlot$id_match+toPlot$id_mc3)
toPlot$perc_MC3_in_PCAWG = toPlot$id_match/(toPlot$id_match+toPlot$id_pcawg)

#Add match color of 25char TCGA barcode match 
iddat <- fread(file=args[2],header=TRUE,sep=" ")
iddat$ord = rowSums(iddat[,3:7])

iddat$char12 = substr(iddat$EXOME,1,12)
toPlot$char12 = substr(toPlot$mc3_exome_barcode,1,12)


toPlot2 <- merge(toPlot,iddat,by="char12")




p <- ggplot(toPlot2, aes(x=perc_PCAWG_in_MC3, y=perc_MC3_in_PCAWG, color=factor(ord)))
p <- p + geom_point(shape=16,alpha=.7)
p <- p + scale_x_continuous(labels = percent)
p <- p + scale_y_continuous(labels = percent)
p <- p + guides(color=guide_legend(reverse=TRUE))
p <- p + theme_bw()
p <- p + theme(panel.border = element_blank())
p <- p + theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))
p <- p + theme(legend.position=c(0,1), legend.justification=c(0,1))
p <- p + xlab("Matched/(Matched+MC3) variants") + ylab("Matched/(Matched+ICGC) variants")
p <- p + geom_vline(xintercept = .8,colour="#ff0066", linetype = "longdash")
p <- p + geom_hline(yintercept = .8,colour="#ff0066", linetype = "longdash")
p <- p + scale_color_manual(values=brewer.pal(name="Set1",n=5)[2:5])
p1 <- ggExtra::ggMarginal(p,type = 'histogram', margins = 'both',size=10,xparams=list(fill="#5A80A6"),yparams=list(fill="#BE312D"))
p1


pdf(args[3],width=4,height=4,useDingbats=F)
print(p1)
dev.off()



q1 <- dim(toPlot[which(toPlot2$perc_PCAWG_in_MC3 >= 0.8 & toPlot$perc_MC3_in_PCAWG >= 0.8),])
q2 <- dim(toPlot[which(toPlot2$perc_PCAWG_in_MC3 < 0.8 & toPlot$perc_MC3_in_PCAWG >= 0.8),])
q3 <- dim(toPlot[which(toPlot2$perc_PCAWG_in_MC3 < 0.8 & toPlot$perc_MC3_in_PCAWG < 0.8),])
q4 <- dim(toPlot[which(toPlot2$perc_PCAWG_in_MC3 >= 0.8 & toPlot$perc_MC3_in_PCAWG < 0.8),])

tot <- dim(toPlot)

q <- rbind(q1,q2,q3,q4,tot)
write(q,args[4],sep="\t")

