library(data.table)
library(likert)
library(plyr)
library(grid)
library(dplyr)

#On my system this requires conda env likert 

#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

dat = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

dat$match = ifelse(!is.na(dat$Chromosome) & !is.na(dat$"Chromosome:1"),2,"NA")
dat$match = ifelse(is.na(dat$Chromosome) & !is.na(dat$"Chromosome:1"),3,dat$match)
dat$match = ifelse(is.na(dat$"Chromosome:1") & !is.na(dat$Chromosome),1,dat$match)


a <- data.frame("FilterName"="Together", "Value"=dat$match)

#This is MUSE and NOTMUSE
muse <- dat[which(grepl("MUSE", dat$CENTERS) | grepl("muse",dat$i_Callers)),]
b <- data.frame("FilterName"="MUSE","Value"=muse$match)
notmuse <- dat[which(!grepl("MUSE", dat$CENTERS) & !grepl("muse",dat$i_Callers)),]
c <- data.frame("FilterName"="notMUSE","Value"=notmuse$match)

#This is BROAD and MUTECT and not these 
mutect <- dat[which(grepl("MUTECT", dat$CENTERS) | grepl("broad",dat$i_Callers)),]
e <- data.frame("FilterName"="MUTECT_BROAD","Value"=mutect$match)
notmutect <- dat[which(!grepl("MUTECT", dat$CENTERS) & !grepl("broad",dat$i_Callers)),]
f <- data.frame("FilterName"="notMUTECT_BROAD","Value"=notmutect$match)

#Variants not found by Muse or Mutect 
neith <- dat[which(grepl("MUSE|MUTECT", dat$CENTERS) | grepl("muse|broad",dat$i_Callers)),]
g <- data.frame("FilterName"="MUSE_MUTECT","Value"=neith$match)
notneith <- dat[which(!grepl("MUSE|MUTECT", dat$CENTERS) & !grepl("muse|broad",dat$i_Callers)),]
h <- data.frame("FilterName"="notMUSE_MUTECT","Value"=notneith$match)


d <- rbind(a,b,c,e,f,g,h)
d$ID = seq(1:dim(d)[1])
dat2 <- dcast(d,d$ID ~ d$FilterName, value.var="Value")
dat2$"d$ID" <- NULL

l <- likert(dat2,nlevels=3)
pdf(args[2],height=3,width=7)
plot(l,include.histogram=TRUE,plot.missing=F)
dev.off()
