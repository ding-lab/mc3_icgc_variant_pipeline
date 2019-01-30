library(data.table)
library(reshape2)
library(sjPlot)
library(ggplot2)
library(RColorBrewer)



args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


dat <- fread(args[1],header=F)

colnames(dat) <- c("FilterName","Value")
dat$ID = seq(1:dim(dat)[1])

dat2 <- dcast(dat,dat$ID ~ dat$FilterName, value.var="Value")
dat2$"dat$ID" <- NULL
#dat2$"nonpreferredpair" <- NULL 
sjp.likert(dat2, values = "sum.outside",sort.frq = "neg.asc",show.prc.sign = TRUE,geom.colors = brewer.pal(3,"Set1"),)


pdf(args[2],height=10,width=7,,useDingbats=F)
sjp.likert(dat2, values = "sum.outside",sort.frq = "neg.asc",show.prc.sign = TRUE,geom.colors = brewer.pal(3,"Set1"),)
dev.off()


