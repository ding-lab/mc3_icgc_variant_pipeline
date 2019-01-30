library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(scales)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

dat = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA")
#colnames(dat)

#table(data.frame(dat$ICGC_Cluster1,dat$clonal.ix))

ui <- dat[which(is.na(dat$Chromosome)),]
#table(ui$ICGC_Cluster1)


um <- dat[which(is.na(dat$"Chromosome:1")),]
#table(um$clonal.ix)

match <- dat[which(!is.na(dat$Chromosome) & !is.na(dat$"Chromosome:1") ),]


p <- ggplot(match,aes(y=i_VAF,x=ICGC_Cluster1))
p <- p+geom_boxplot()
p

p <- ggplot(match,aes(y=t_alt_count/t_depth,x=clonal.ix))
p <- p+geom_boxplot()
p

match$clusternum <- as.numeric((apply(match[,176:180],1,which.max)))
dat$clusternum <- as.numeric((apply(dat[,176:180],1,which.max)))

p <- ggplot(match,aes(y=t_alt_count/t_depth,x=factor(clusternum)))
p <- p+geom_boxplot()
p

p <- ggplot(match,aes(y=i_VAF,x=factor(clusternum)))
p <- p+geom_boxplot()
p


#### MAKE THIS PLOT FOR ICGC ####

ICGC_dat <- dat[which(!is.na(dat$"Chromosome:1")),]
ICGC_dat$ICGC_Uniq = ifelse(is.na(ICGC_dat$Chromosome),"Unique","Matched")
dat2 <- ICGC_dat %>%
  group_by(clusternum) %>% mutate(n = n()) %>%
  mutate(label = paste0(clusternum,'\nN = ',n))

dat2$n <- NULL

summary = dat2 %>% group_by(label,ICGC_Uniq) %>%
  tally %>%
  group_by(label) %>%
  mutate(pct = n/sum(n),
         n.pos = cumsum(n) - 0.5)

#TODO there are 4178 mutations that are have unassigned subclonal structure, mark this in the file and we will have to go back to these numbers to adjust it, one example is a SNV frame_shift_mutations for  301d6ce3-4099-4c1d-8e50-c04b7ce91450 at position 3:156396062 according to the MAF, but the other has a mutations at 3:156396061 in the indel section of the cluster assignments file. Make sure to not this in the figure legend. This may extend further into that mutation clonality match. 

p <- ggplot(summary[which(summary$label != "NA\nN = 3958"),], aes(x=label, y=pct, fill=ICGC_Uniq))
p <- p+geom_bar(stat="identity")
p <- p+geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%"), y=ifelse(pct>.5,.5,pct*.5)),colour="white")
p <- p+labs(y="Percent of mutations",x="Subclonal designation")
p <- p+scale_y_continuous(labels=percent)
p <- p+scale_fill_manual(values=brewer.pal(name="Set2",n=3))
p <- p+theme_minimal()
p <- p+theme(legend.position="bottom")
p


pdf(args[2],height=4,width=4,useDingbats=F)
print(p)
dev.off()

ICGC_summary <- summary



#TODO, and I'm talking to eduard about this one: There are 59,232 mutations without TCGA clonality extimations in these data. Most of these come from a handful of samples but still these were reduced, but hopefully I will be able to recover this for a down stream anlayis., Keep in mind these mutations were restricted to SNV called PCAWG variants. 

#### Now I need to make a similar plot for TCGA/MC3
TCGA_dat <- dat[which(!is.na(dat$"Chromosome")),]
TCGA_dat$TCGA_Uniq = ifelse(is.na(TCGA_dat$"Chromosome:1"),"Unique","Matched")
dat2 <- TCGA_dat %>%
  group_by(clonal.ix) %>% mutate(n = n()) %>%
  mutate(label = paste0(clonal.ix,'\nN = ',n))

dat2$n <- NULL

summary = dat2 %>% group_by(label,TCGA_Uniq) %>%
  tally %>%
  group_by(label) %>%
  mutate(pct = n/sum(n),
         n.pos = cumsum(n) - 0.5)


#Here I just want to change the words so that they appear in order 
summary$label = ifelse(summary$label == "FALSE\nN = 26406", "Sub-clonal\nN = 26406",summary$label)
summary$label = ifelse(summary$label == "TRUE\nN = 87433", "Clonal\nN = 87433",summary$label)


p <- ggplot(summary[which(summary$label != "NA\nN = 58669"),], aes(x=label, y=pct, fill=TCGA_Uniq))
# <- ggplot(summary, aes(x=label, y=pct, fill=TCGA_Uniq))
p <- p+geom_bar(stat="identity")
p <- p+geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%"), y=ifelse(pct>.5,.5,pct*.5)),colour="white")
p <- p+labs(y="Percent of mutations",x="Subclonal designation")
p <- p+scale_y_continuous(labels=percent)
p <- p+scale_fill_manual(values=brewer.pal(name="Set2",n=3))
p <- p+theme_minimal()
p <- p+theme(legend.position="bottom")
p

pdf(args[3],height=4,width=4,useDingbats=F)
print(p)
dev.off()


