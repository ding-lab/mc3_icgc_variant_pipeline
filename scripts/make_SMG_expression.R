library(data.table)
library(ggplot2)
library(doParallel)


dat <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv")
#can <- fread("FUT9.LIHC.cis.txt",header=F)
#can <- fread("FUT9.LUSC.cis.txt",header=F)
#can <- fread("SFTPB.LUAD.cis.txt",header=F) #SIGNIFICANT
#can <- fread("SNHG14.LUSC.cis.txt",header=F) 
#can <- fread("MMP16.LUSC.cis.txt",header=F) 
#can <- fread("AFF2.STAD.cis.txt",header=F) 
#can <- fread("PGR.SKCM.cis.txt",header=F)
#can <- fread("ERBB4.UCEC.cis.txt",header=F)


etest <- function(x,y){
    gene <- dat[which(dat$gene_id == x),]
    can <- fread(y,header=F)
    tgene <- data.frame(t(gene))
    tgene$char12 = substr(rownames(tgene),1,12)
    tg2 <- tgene[-1,]
    can$char12 <- substr(can$V1,1,12)
    totest <- merge(tg2,can,by="char12")
    return(totest)
}

allgenes <- unique(dat$gene_id)


cl <- makeCluster(60) #Number
registerDoParallel(cl)

global_t <- foreach( i=1:length(allgenes), .combine=rbind, .export = 'fread') %dopar%{
    x = allgenes[i]
    y = "ERBB4.UCEC.cis.txt"
    look <- etest(x,y)
    a <- look[which(look$V2 == "MUT"),]
    b <- look[which(look$V2 == "WT"),]
    if(dim(a)[1] > 3 & dim(b)[1]>3 & mean(as.numeric(as.character(a$t.gene.))) != mean(as.numeric(as.character(b$t.gene.)))) {
        d <- t.test(as.numeric(as.character(a$t.gene.)),as.numeric(as.character(b$t.gene.)))
        e <- data.frame("Gene"=i,"PVAL"=d$p.value)
        return(e)
    }
}
stopCluster(cl)

global_t$Hugo = allgenes[global_t$Gene]

head(global_t[order(global_t$PVAL),],20)


write.table(global_t,'bestERBB4.UCEC.optim.txt',quote=F,sep="\t",row.names=F)



#################So I looked at a couple non-coding things and didn't find

#### BEST OF  SFTPB #####

sftpb <- fread("bestSFTPB.LUAD.optim.v2.txt")
head(sftpb[order(sftpb$PVAL),],20)


x = "EGFL6|25975"
y = "SFTPB.LUAD.cis.txt"

look <- etest(x,y)

p <- ggplot(look,aes(x=V2,y=as.numeric(as.character(t.gene.))))
p <- p + geom_dotplot(binaxis="y",stackdir="center")
p


