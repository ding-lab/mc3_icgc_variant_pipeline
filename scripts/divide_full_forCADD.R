library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

full = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))

mc3 <- unique(data.frame("CHR"=full$Chromosome, "POS"=full$Start_Position, "END"=full$End_Position, "VAR"=paste(full$Reference_Allele,"/",full$Tumor_Seq_Allele2,sep=""), "STRAND"= full$Strand, "Varclass"=full$Variant_Type))
mc3snps <- mc3[which(mc3$Varclass == "SNP"),]

pcawg <- unique(data.frame("CHR"=full$"Chromosome:1", "POS"=full$"Start_position:1", "END"=full$"End_position:1", "VAR"=paste(full$"Reference_Allele:1","/",full$"Tumor_Seq_Allele2:1",sep=""), STRAND=full$"Strand:1", "Varclass"=full$"Variant_Type:1"))
pcawgsnps <- pcawg[which(pcawg$Varclass == "SNP"),]
psudoVCFpcawg <- pcawgsnps[,1:5]
psudoVCFmc3 <- mc3snps[,1:5]
psudoVCF <- unique(rbind(psudoVCFpcawg,psudoVCFmc3))
write.table(psudoVCF,args[5],sep="\t",row.names=F,quote=F)

ps <- psudoVCF
chrOrder <-c((1:22),"X","Y","M")
ps$CHR <- factor(ps$CHR, chrOrder, ordered=TRUE)
sps <- ps[order(ps$CHR,ps$POS),]
ssps <- split(sps, rep(1:3, length.out = nrow(sps), each = ceiling(nrow(sps)/3)))
write.table(ssps[[1]],args[2],sep="\t",quote=F,row.names=F,col.names=FALSE)
write.table(ssps[[2]],args[3],sep="\t",quote=F,row.names=F,col.names=FALSE)
write.table(ssps[[3]],args[4],sep="\t",quote=F,row.names=F,col.names=FALSE)


