library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}


#################################
#CODE written by William Meyerson 
#Yale University
#################################


options(stringsAsFactors = F)
x = as.data.frame(fread(args[1]))


#Calculate Total reads 
x$t_tot = x$t_alt_count + x$t_ref_count
x$p_tot = x$`t_alt_count:1` + x$`t_ref_count:1`
#Calculate VAF
x$t_VAF = x$t_alt_count/x$t_tot
x$p_VAF = x$`t_alt_count:1`/x$p_tot
#Find unique and matched calls 
x$t_only = !is.na(x$Tumor_Sample_Barcode) & is.na(x$Project_Code)
x$p_only = is.na(x$Tumor_Sample_Barcode) & !is.na(x$Project_Code)
x$both = !is.na(x$Tumor_Sample_Barcode) & !is.na(x$Project_Code)

#Rounding-up Mechanims
x$t_VAFbin = ceiling(x$t_VAF*50)/50
x$p_VAFbin = ceiling(x$p_VAF*50)/50
#X is 203125 rows at this point

x = x[unique(c(which(x$t_alt_count>=3), which(x$`t_alt_count:1`>=3))),]
#x is 202755 rows at this point less than or equal to 3

# Have to filter out variants with altcount < 3
x$both[which(x$t_alt_count < 3)] = F
x$both[which(x$`t_alt_count:1` < 3)] = F

#Make sure all read counts are finite
x$m_tot = apply(x[,c("p_tot", "t_tot")], 1, max, na.rm=T)
x = x[!is.infinite(x$m_tot),]
#x is 202751 right now. 

#This make the transition transversion chars
x$varO = paste(x$Tumor_Seq_Allele1, x$Tumor_Seq_Allele2, sep=">")
x$varO[grep("NA", x$varO)] = paste(x[grep("NA", x$varO),"Tumor_Seq_Allele1:1"], x[grep("NA", x$varO), "Tumor_Seq_Allele2:1"], sep=">")

goodVarO = c("A>C", "A>G", "A>T", "C>G", "C>T", "G>A", "G>C", "T>A", "T>C","T>G")

x = x[which(x$varO %in% goodVarO),]
#x is now 156703 rows

# q = quantile(x$m_tot)

# Tukey's rule says that the outliers are values more than 1.5 times the interquartile range from the quartiles

# iq = as.numeric(q[4] - q[2])
# tukhigh = q[4] + 1.5*iq
# 
# x = x[which(x$m_tot <= tukhigh),]

#This bit of code squishes VAF bi
ag.vaf.emp = aggregate(both ~ t_VAFbin, FUN=mean, data=x)
ag.vaf.pemp = aggregate(both ~ p_VAFbin, FUN=mean, data=x)
ag.cnt.pemp = aggregate(Tumor_Sample_Barcode ~ p_VAFbin, FUN=length, data=x)

y = as.data.frame(fread(input=args[2], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139))))
y$t_tot = y$t_alt_count + y$t_ref_count
y$p_tot = y$`t_alt_count:1` + y$`t_ref_count:1`
y$t_bin1MB = paste(y$Chromosome, ceiling(y$Start_Position/10^6), sep="_")
y$p_bin1MB = paste(y$`Chromosome:1`, ceiling(y$`Start_position:1`/10^6), sep="_")
y.ag.t_tot.t_bin1MB = aggregate(t_tot ~ t_bin1MB, FUN=quantile, 0.5, data=y)
y.ag.p_tot.p_bin1MB = aggregate(p_tot ~ p_bin1MB, FUN=quantile, 0.5, data=y)
names(y.ag.t_tot.t_bin1MB)[1] = "bin1MB"
names(y.ag.p_tot.p_bin1MB)[1] = "bin1MB"
y.ag.m_tot.m_bin1MB = merge(y.ag.t_tot.t_bin1MB, y.ag.p_tot.p_bin1MB, all=T)
y.ag.m_tot.m_bin1MB[is.na(y.ag.m_tot.m_bin1MB)] = 0
# plot(y.ag.m_tot.m_bin1MB$t_tot, y.ag.m_tot.m_bin1MB$p_tot, xlab="Median covergage TCGA",
     # ylab="Median coverage PCAWG")

# splitY = split(y, y$mc3_exome_barcode)


####


# simulate VAFs with coverage
x$t_bin1MB = paste(x$Chromosome, ceiling(x$Start_Position/10^6), sep="_")
x$p_bin1MB = paste(x$`Chromosome:1`, ceiling(x$`Start_position:1`/10^6), sep="_")
x$m_bin1MB = x$t_bin1MB
x$m_bin1MB[grep("NA", x$t_bin1MB)] = x$p_bin1MB[grep("NA", x$t_bin1MB)]

u.pbin = y.ag.m_tot.m_bin1MB$bin1MB

x$p_simReads = NA
for(i in 1:length(u.pbin)) {
  myBin = u.pbin[i]
  ix = which(x$m_bin1MB == myBin)
  pot_reads = x$p_tot
  clear_reads = pot_reads[!is.na(pot_reads)]
  if(length(clear_reads)==0) {
    clear_reads = 0
  }
  x$p_simReads[ix] = sample(clear_reads, length(ix), replace=T)
}


wdf = x[which(x$t_alt_count >=3),c("t_alt_count", "t_tot", "p_simReads")]
wdf$t_VAF_laplace = (wdf$t_alt_count+1)/(wdf$t_tot + 2)


p.seq = rep(NA, nrow(wdf))
for(i in 1:nrow(wdf)) {
  p.seq[i] = rbinom(n=1, size=wdf$p_simReads[i], prob=wdf$t_VAF_laplace[i])
}
wdf$p.mock = p.seq >= 3
wdf$t_VAFbin = ceiling(wdf$t_alt_count/wdf$t_tot*50)/50
ag.vaf.sim = aggregate(p.mock ~ t_VAFbin, FUN=mean, data=wdf)


# plot(ag.vaf.emp$t_VAFbin, ag.vaf.emp$both, ylim=c(0,1), pch=16)
# points(ag.vaf.pemp$p_VAFbin, ag.vaf.pemp$both, col="blue", pch=16)
# points(ag.vaf.sim$t_VAFbin, ag.vaf.sim$p.mock, col="purple", pch=16)




######Simulated it out to another dataset 

x$t_simReads = NA
for(i in 1:length(u.pbin)) {
  myBin = u.pbin[i]
  ix = which(x$m_bin1MB == myBin)
  pot_reads = x$t_tot
  clear_reads = pot_reads[!is.na(pot_reads)]
  if(length(clear_reads)==0) {
    clear_reads = 0
  }
  x$t_simReads[ix] = sample(clear_reads, length(ix), replace=T)
}


tdf = x[which(x$`t_alt_count:1` >=3 & !is.na(x$p_tot)),c("t_alt_count:1", "p_tot", "t_simReads")]
tdf$p_VAF_laplace = (tdf$`t_alt_count:1`+1)/(tdf$p_tot + 2)

t.seq = rep(NA, nrow(tdf))
for(i in 1:nrow(tdf)) {
  t.seq[i] = rbinom(n=1, size=tdf$t_simReads[i], prob=tdf$p_VAF_laplace[i])
}
tdf$t.mock = t.seq >= 3
tdf$p_VAFbin = ceiling(tdf$`t_alt_count:1`/tdf$p_tot*50)/50
ag.vaf.tsim = aggregate(t.mock ~ p_VAFbin, FUN=mean, data=tdf)


pdf(args[3], useDingbats = F)
plot(ag.vaf.emp$t_VAFbin, ag.vaf.emp$both, ylim=c(0,1), pch=16,
     xlab="VAF in consortium X", ylab="Recovery fraction in consortium Y")
points(ag.vaf.pemp$p_VAFbin+0.005, ag.vaf.pemp$both, col="blue", pch=16, cex=1)
points(ag.vaf.sim$t_VAFbin, ag.vaf.sim$p.mock, col="orange", pch=16)
points(ag.vaf.tsim$p_VAFbin+0.005, ag.vaf.tsim$t.mock, col="green", pch=16, cex=1)
legend(0.4, 0.6, legend=c("Simulated, PCAWG = X", "Simulated, TCGA = X",
                          "Observed, PCAWG = X", "Observed, TCGA = X"),
       col=c("green", "orange", "black", "blue"), lty=1, cex=1)
dev.off()

