library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)



args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

df_variant = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

gene_model_pth <- args[2]
df_raw <- fread(gene_model_pth)

# Merge the two
df_merged = cbind(df_variant[1:nrow(df_raw), ], df_raw)


# Lets take this a little bit of a different direction: 
df_merged$match = ifelse(!is.na(df_merged$Chromosome) & !is.na(df_merged$"Chromosome:1"),1,0)
df_merged$PCAWG_only = ifelse(is.na(df_merged$Chromosome) & !is.na(df_merged$"Chromosome:1"),1,0)
df_merged$MC3_only = ifelse(is.na(df_merged$"Chromosome:1") & !is.na(df_merged$Chromosome),1,0)


uu_p <- fread(args[3])

uu_p$line <- do.call(paste,as.list(uu_p[,111:124]))
df_merged$line <- do.call(paste,as.list(df_merged[,111:124]))

uu_pcawg <- df_merged[which(df_merged$line %in% uu_p$line),]

#uu_mc3 <- df_merged[which(df_merged$MC3_only == 1 & (df_merged$FILTER == "PASS" | df_merged$FILTER == "wga")),]
#uu_match <- df_merged[which(df_merged$match == 1 & (df_merged$FILTER == "PASS" | df_merged$FILTER == "wga" | is.na(df_merged$FILTER))),]

uu_mc3 <- df_merged[which(df_merged$MC3_only == 1),]
uu_match <- df_merged[which(df_merged$match == 1 ),]


uu_pcawg$UNIQUE = "PCAWG"
uu_mc3$UNIQUE = "MC3"
uu_match$UNIQUE = "MATCH"

all <- rbind(uu_mc3,uu_pcawg,uu_match)
all$line <- NULL

look <- data.frame(all %>% group_by(Annotation,UNIQUE,mycut = cut(all$'Normalized position', breaks = 15)) %>% tally())

look$annoOrder = factor(look$Annotation,level=c("5'UTR","exon1","exon2","exon3","3'UTR"))
look2 <- look[which(!is.na(look$Annotation)),]




#NOTE This has room for filtering: 
#By proportion: 
gmp <- ggplot(look2,aes(x=mycut,y=n,fill=UNIQUE))
gmp <- gmp+geom_bar(stat="identity",position="fill")
gmp <- gmp+facet_grid(~annoOrder )
gmp <- gmp+theme_minimal()
gmp <- gmp+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
gmp

#By count
gmn <- ggplot(look2,aes(x=mycut,y=n,fill=UNIQUE))
gmn <- gmn+geom_bar(stat="identity")
gmn <- gmn+facet_grid(~annoOrder )
gmn <- gmn+theme_minimal()
gmn <- gmn+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
gmn


pdf(args[4],height=4,width=9,useDingbats=F)
print(gmp)
dev.off()


pdf(args[5],height=4,width=9,useDingbats=F)
print(gmn)
dev.off()


