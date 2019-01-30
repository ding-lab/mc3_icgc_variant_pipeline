library(data.table)
library(dplyr)
library(ggplot2)



#For testing 
#args = c("processed_data/removed.by.coverage.genome.maf","/diskmnt/Projects/ICGC_MC3/Reduction_Beds/pseudogenes.txt","/diskmnt/Projects/ICGC_MC3/Data/cancer.gene.299.txt","figures/UTR3.Cancer.Histogram.pdf","figures/UTR5.Cancer.Histogram.pdf","figures/MISS.Cancer.Histogram.pdf","figures/Coverage_sunburst.pdf")
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


dat <- fread(args[1],header=F,na.strings="NA",colClasses=list(character=c(17,33,35)))
pseudo <- fread(args[2],header=F)

'%!in%' <- function(x,y)!('%in%'(x,y))

ps <- (dat[which(dat$V43 %in% pseudo$V1)])# = 2890
dat2 <- (dat[which(dat$V43 %!in% pseudo$V1)])
typ <- data.frame(dat2 %>% group_by(V5) %>% tally())


#By cancer genes 299 
can <- fread(args[3],header=F)
dat3 <- dat2[which(dat2$V43 %in% can$V1),]

genic <- data.frame(dat3 %>% group_by(V43,V5) %>% tally())



#FOR 3'UTRs 
UTR3 <- genic[which(genic$V5 == "3'UTR"),]
Top15_UTR <- head(UTR3[order(UTR3$n,decreasing=T),],15)

p <- ggplot(Top15_UTR, aes(x=reorder(V43,-n),y=n))
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p <- p + ylab("Variant count")
p <- p + xlab("Gene")
p <- p+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=.5))
p

pdf(args[4],height=3, width=4,useDingbats=F)
print(p)
dev.off()

#For 5'Flank
UTR5 <- genic[which(genic$V5 == "5'UTR"),]
Top15_UTR5 <- head(UTR5[order(UTR5$n,decreasing=T),],15)

p <- ggplot(Top15_UTR5, aes(x=reorder(V43,-n),y=n))
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p <- p + ylab("Variant count")
p <- p + xlab("Gene")
p <- p+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=.5))
p

pdf(args[5],height=3, width=4,useDingbats=F)
print(p)
dev.off()


#Missense
MISS = genic[which(genic$V5 == "Missense_Mutation"),]
Top15_MISS <- head(MISS[order(MISS$n,decreasing=T),],15)

p <- ggplot(Top15_MISS, aes(x=reorder(V43,-n),y=n))
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p <- p + ylab("Variant count")
p <- p + xlab("Gene")
p <- p+theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=.5))
p

pdf(args[6],height=3, width=4,useDingbats=F)
print(p)
dev.off()



##### THE ACTUAL SUNBURST ###### 
#### THIS STILL NEEDS TO BE CODED UP###### 

lev2 = data.frame(dat2 %>% group_by(V5) %>% tally())

lev3 = data.frame(dat3 %>% group_by(V5) %>% tally()) 

covered = 387166-dim(dat)[1]

FLANK = c("3'UTR","5'Flank","5'UTR")
flank_all = sum(lev2[which(grepl("'",lev2$V5)),]$n)
flank_can = sum(lev3[which(grepl("'",lev3$V5)),]$n)
flank_lev2 = flank_all - flank_can

INTRON_IGR = c("IGR","Intron")
intronIGR_all = sum(lev2[which(lev2$V5 %in% INTRON_IGR),]$n)
intronIGR_can = sum(lev3[which(lev3$V5 %in% INTRON_IGR),]$n)
intronIGR_lev2 = intronIGR_all-intronIGR_can

RNA_lncRNA = c("lincRNA","RNA")
rnaLNC_all = sum(lev2[which(lev2$V5 %in% RNA_lncRNA),]$n)
#rnaLNC_can = sum(lev3[which(lev3$V5 %in% RNA_lncRNA),]$n)
#rnaLNC_lev2 = rnaLNC_all-rnaLNC_can

GENE = c("Frame_Shift_Ins","Frame_Shift_Del","Splice_Site","Silent","Missense_Mutation","In_Frame_Del","Nonsense_Mutation","De_novo_Start_InFrame","De_novo_Start_OutOfFrame","Start_Codon_SNP","Start_Codon_Del","In_Frame_Ins","Start_Codon_Del","Stop_Codon_Del")
gene_all = sum(lev2[which(lev2$V5 %in% GENE),]$n)
gene_can = sum(lev3[which(lev3$V5 %in% GENE),]$n)
gene_lev2 = gene_all-gene_can


pseudo = dim(ps)[1]



values = c(covered,flank_lev2,flank_can,intronIGR_lev2,intronIGR_can,rnaLNC_all,gene_lev2,gene_can,pseudo)



df <- data.frame(
    'level1' = c('Covered', 'Not covered', 'Not covered', 'Not covered', 'Not covered', 'Not covered', 'Not covered', 'Not covered', 'Not covered'),
    'level2' = c('Covered', 'Flank', 'Flank', 'Introns & IGR', 'Introns & IGR', 'RNA & lncRNA', 'Genes', 'Genes', 'Pseudo'),
    'level3' = c('Covered_can', 'Flank_data', 'Flank_can', 'Introns & IGRdata', 'Introns & IGRcan', 'RNA & lncRNA_data', 'Genes_data', 'Genes_can', 'Pseudo'),
    'value' = values)

p <- ggplot(df, aes(y=value)) +
    geom_bar(aes(fill=level1, x=0), width=.5, stat='identity') +
    geom_bar(aes(fill=level2, x=.25), width=.25, stat='identity') +
    geom_bar(aes(fill=level3, x=.5), width=.25,stat="identity") +
    coord_polar(theta='y') +
    theme_classic()


pdf(args[7], height=8, width=8, useDingbats = F)
print(p)
dev.off()



