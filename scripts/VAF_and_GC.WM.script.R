##########
## CODE authored by William Meyerson of Gerstein Lab at Yale University
## Contact william.meyerson@yale.edu [until 05/2020] or william.ulysses@gmail.com [thereafter] with questions
##########

##### Table of Contents
## 1. Set up environment
## 2. Load and prepare inputs
## 3. Idealized VAF model
## 4. Empirical VAF model
## 5. GC content

### 1. Set up environment

options(stringsAsFactors = F)
list.of.packages <- c("data.table", "ranger", "broom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("data.table")
library("ranger")
library("broom")

### 2. Load and prepare inputs
# Obtain file paths
# Load in files
# Remove variants that fail chosen filters
# For certain variables, extract whichever of PCAWG or MC3 is non-NA
# Extract only SNVs
# Obtain patient id
# Exclude hypermutator samples
# Standardize read depths and variant allele frequencies
# Designate variant capture status
# Require variants to have at least three alternate reads in the tumor

# Obtain file paths
MAF = commandArgs(trailingOnly = T)[1]
PURITY.TCGA = commandArgs(trailingOnly = T)[2]
PURITY.PCAWG = commandArgs(trailingOnly = T)[3]
CADD = commandArgs(trailingOnly = T)[4]

###@@@### Development only ###@@@###
MAF = "~/Downloads/full6.tsv"
PURITY.PCAWG = "/Users/ulysses/Documents/R/VALID/purity_ploidy.txt"
PURITY.TCGA = "/Users/ulysses/Downloads/tcgaPurity.csv"
CADD = "~/Documents/R/bailCad/outBailCad.txt"
###@@@########################@@@###

# Load in files
maf = as.data.frame(fread(MAF))
purity = read.csv(PURITY.PCAWG, sep="\t")[,1:2]
tcgaPurity = read.csv(PURITY.TCGA, skip=3)
cadd = as.data.frame(fread(CADD))
caseKey = read.csv(CASE) 

# Remove variants that fail chosen filters
maf = maf[-unique(c(grep("oxog", maf$FILTER),
                grep("nonpreferredpair", maf$FILTER))),]


# For certain variables, extract whichever of PCAWG or MC3 is non-NA
pickFinite = function(S) {
  s2 = S[which(!is.na(S))]
  if(length(s2) > 1) {
    s2 = sample(s2, 1)
  }
  return(s2)
}

maf$my.chr = apply(maf[,c("Chromosome", "Chromosome:1")], 1, pickFinite)
maf$my.pos = apply(maf[,c("Start_Position", "Start_position:1")], 1, pickFinite)
maf$my.ref = apply(maf[,c("Tumor_Seq_Allele1", "Tumor_Seq_Allele1:1")], 1, pickFinite)
maf$my.alt = apply(maf[,c("Tumor_Seq_Allele2", "Tumor_Seq_Allele2:1")], 1, pickFinite)

# Extract only SNVs
ACGT = c("A", "C", "G", "T")
maf = maf[which(maf$my.ref %in% ACGT & maf$my.alt %in% ACGT),]

# Obtain patient id
maf$first16 = gsub("(....-..-....-...).*", "\\1", maf$mc3_exome_barcode)

# Exclude hypermutator samples
umc3 = unique(maf$mc3_exome_barcode)
hypermutators = c("TCGA-AA-3977-01A-01W-0995-10", "TCGA-AA-A00N-01A-02D-A17O-10", "TCGA-AD-6964-01A-11D-1924-10", "TCGA-BR-6452-01A-12D-1800-08", 
           "TCGA-CA-6717-01A-11D-1835-10", "TCGA-CA-6718-01A-11D-1835-10", "TCGA-F5-6814-01A-31D-1924-10", "TCGA-GN-A266-06A-11D-A197-08")
maf = maf[!maf$mc3_exome_barcode %in% hypermutators,]


# Standardize read depths and variant allele frequencies
maf$t_tot = maf$t_alt_count + maf$t_ref_count
maf$p_tot = maf$`t_alt_count:1` + maf$`t_ref_count:1`
maf$t_VAF = maf$t_alt_count/maf$t_tot
maf$p_VAF = maf$`t_alt_count:1`/maf$p_tot
maf$m_tot = apply(maf[,c("p_tot", "t_tot")], 1, max, na.rm=T)
maf = maf[!is.infinite(maf$m_tot),]

# Designate variant capture status
maf$t_only = maf$MC3_only
maf$p_only = maf$PCAWG_only
maf$both = maf$match


# Require variants to have at least three alternate reads in the tumor [turn off]
# maf = maf[unique(c(which(maf$t_alt_count>=3), which(maf$`t_alt_count:1`>=3))),]
# maf$both[which(maf$t_alt_count < 3)] = F
# maf$both[which(maf$`t_alt_count:1` < 3)] = F


### 3. Idealized VAF model

### Need wdf, tdf, maf, df.ubinz, ag.vaf.emp, ag.vaf.pemp, ag.vaf.sim, ag.vaf.tsim


# Assign variants to VAF bins
THR = 50 # How many VAF bins?
maf$t_VAFbin = ceiling(maf$t_VAF*THR)/THR
maf$p_VAFbin = ceiling(maf$p_VAF*THR)/THR

# Tablulate co-call rate of variants by MC3 VAF bin and PCAWG VAF bin
ag.vaf.emp = aggregate(both ~ t_VAFbin, FUN=mean, data=maf)
ag.vaf.pemp = aggregate(both ~ p_VAFbin, FUN=mean, data=maf)
ag.cnt.pemp = aggregate(Tumor_Sample_Barcode ~ p_VAFbin, FUN=length, data=maf)

upt = unique(maf$mc3_exome_barcode)

## Simulate the PCAWG read depth for each variant in TCGA
## By sampling from the PCAWG read depths of other variants called in the involved patient
maf$p_simReads = NA
for(i in 1:length(upt)) {
  myPt = upt[i]
  ix = which(maf$mc3_exome_barcode == myPt)
  pot_reads = maf$p_tot[ix]
  clear_reads = pot_reads[!is.na(pot_reads)]
  if(length(clear_reads)==0) {
    clear_reads = 0
  }
  maf$p_simReads[ix] = sample(clear_reads, length(ix), replace=T)
}
maf$p_simReads[!is.na(maf$p_tot)] = maf$p_tot[!is.na(maf$p_tot)]

wdf = maf[which(maf$t_alt_count >=3),c("t_alt_count", "t_tot", "p_simReads", "first16")]


# Use laplace's law of succession to smooth variant allele frequency estimates when alt_counts are extremely high or low
wdf$t_VAF_laplace = (wdf$t_alt_count+1)/(wdf$t_tot + 2)

# Adjust variant allele frequency estimates estimates for sample purity
us = unique(maf[,c("first16", "Tumor_Sample_Barcode:1")])
us = us[complete.cases(us),]
wdf$TSB = us$`Tumor_Sample_Barcode:1`[match(wdf$first16, us$first16)]
wdf$purityTcga = tcgaPurity$ESTIMATE[match(wdf$first16, tcgaPurity$Sample.ID)]
wdf$purityPcawg = purity$purity[match(wdf$TSB, purity$sample)]

wdf$t_VAF.purityAdjusted = wdf$t_VAF_laplace/wdf$purityTcga*wdf$purityPcawg
wdf$t_VAF.purityAdjusted = apply(cbind(wdf$t_VAF.purityAdjusted,1),1,min)
b = maf[!is.na(maf$t_VAF) & !is.na(maf$p_VAF),]
bvaf = merge(aggregate(t_VAF ~ first16, FUN=mean, data=b),
             aggregate(p_VAF ~ first16, FUN=mean, data=b))
bvaf$PtoTrat = bvaf$p_VAF/bvaf$t_VAF
wdf$PtoTrat = bvaf$PtoTrat[match(wdf$first16, bvaf$first16)]
wdf$p_VAF.purityAdjusted = wdf$t_VAF_laplace * wdf$PtoTrat
wdf$p_VAF.purityAdjusted = apply(cbind(wdf$p_VAF.purityAdjusted,1),1,min)
wdf = wdf[!is.na(wdf$p_VAF.purityAdjusted),]

# Simulate alternative allele counts in PCAWG, from estimated PCAWG VAF and estimated PCAWG depth
p.seq = rep(NA, nrow(wdf))
for(i in 1:nrow(wdf)) {
  p.seq[i] = rbinom(n=1, size=wdf$p_simReads[i], prob=wdf$p_VAF.purityAdjusted[i])
}

# Simulate calling variants in simulated PCAWG as variants supported by 3 or more simulated alternate allele counts
wdf$p.mock = p.seq >= 3
wdf$t_VAFbin = ceiling(wdf$t_alt_count/wdf$t_tot*THR)/THR
ag.vaf.sim = aggregate(p.mock ~ t_VAFbin, FUN=mean, data=wdf)

# Repeat for simulating the TCGA read depth, ture VAFs, alt allele counts, and variant call status for each variant in PCAWG
maf$t_simReads = NA
for(i in 1:length(upt)) {
  myPt = upt[i]
  ix = which(maf$mc3_exome_barcode == myPt)
  pot_reads = maf$t_tot[ix]
  clear_reads = pot_reads[!is.na(pot_reads)]
  if(length(clear_reads)==0) {
    clear_reads = 0
  }
  maf$t_simReads[ix] = sample(clear_reads, length(ix), replace=T)
}
maf$t_simReads[!is.na(maf$t_tot)] = maf$t_tot[!is.na(maf$t_tot)]
tdf = maf[which(maf$`t_alt_count:1` >=3 & !is.na(maf$p_tot)),c("t_alt_count:1", "p_tot", "t_simReads", "first16")]
tdf$p_VAF_laplace = (tdf$`t_alt_count:1`+1)/(tdf$p_tot + 2)
tdf$PtoTrat = bvaf$PtoTrat[match(tdf$first16, bvaf$first16)]
tdf$t_VAF.purityAdjusted = tdf$p_VAF_laplace * tdf$PtoTrat
tdf$t_VAF.purityAdjusted = apply(cbind(tdf$t_VAF.purityAdjusted,1),1,min)
tdf = tdf[!is.na(tdf$t_VAF.purityAdjusted),]
t.seq = rep(NA, nrow(tdf))
for(i in 1:nrow(tdf)) {
  t.seq[i] = rbinom(n=1, size=tdf$t_simReads[i], prob=tdf$t_VAF.purityAdjusted[i])
}
tdf$t.mock = t.seq >= 3
tdf$p_VAFbin = ceiling(tdf$`t_alt_count:1`/tdf$p_tot*THR)/THR
ag.vaf.tsim = aggregate(t.mock ~ p_VAFbin, FUN=mean, data=tdf)

## Plot the simulated and observed recovery fraction by VAF bin

pdf("Idealized.VAF.RecoveryRate.pdf", useDingbats = F)
plot(ag.vaf.emp$t_VAFbin, ag.vaf.emp$both, ylim=c(0,1), pch=16,
     xlab="VAF in consortium X", ylab="Recovery fraction in consortium Y", col="red")
points(ag.vaf.pemp$p_VAFbin+0.005, ag.vaf.pemp$both, col="blue", pch=16, cex=1)
points(ag.vaf.sim$t_VAFbin, ag.vaf.sim$p.mock, col="pink", pch=16)
points(ag.vaf.tsim$p_VAFbin+0.005, ag.vaf.tsim$t.mock, col="lightblue", pch=16, cex=1)
legend(0.4, 0.6, legend=c("Simulated, PCAWG = Y", "Simulated, TCGA = Y",
                          "Observed, PCAWG = Y", "Observed, TCGA = Y"),
       col=c("pink", "lightblue", "red", "blue"), lty=1, cex=1)
dev.off()


# The fraction of MC3 variants that are also called in an idealized, VAF-based simulation of PCAWG
print(mean(wdf$p.mock))
# The fraction of MC3 variants that are also called in PCAWG
print(sum(maf$both) / sum(maf$t_only | maf$both))

# The fraction of PCAWG variants that are also called in an idealized, VAF-based simulation of MC3
print(mean(tdf$t.mock))
# The fraction of PCAWG variants that are also called in MC3
print(sum(maf$both) / sum(maf$p_only | maf$both))

# Of the variants in MC3 that are missed in PCAWG, what fraction were expected to be missed based on an idealized VAF model?
print((1-mean(wdf$p.mock))/(1-(sum(maf$both)/sum(maf$both | maf$t_only))))
# Of the variants in PCAWG that are missed in MC3, what fraction were expected to be missed based on an idealized VAF model?
print((1-mean(tdf$t.mock))/(1-(sum(maf$both)/sum(maf$both | maf$p_only))))

### 4. Empirical VAF model

maf$t_VAF.25 = maf$t_VAF
maf$t_VAF.75 = maf$t_VAF
maf$t_VAF.9 = maf$t_VAF
maf$t_VAF.1 = maf$t_VAF

maf$p_VAF.25 = maf$p_VAF
maf$p_VAF.75 = maf$p_VAF
maf$p_VAF.9 = maf$p_VAF
maf$p_VAF.1 = maf$p_VAF

msamp = merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(merge(aggregate(t_VAF ~ mc3_exome_barcode, FUN=mean, data=maf),
                                                                          aggregate(t_VAF.25 ~ mc3_exome_barcode, FUN=quantile, 0.25, na.rm=T, data=maf)),
                                                                    aggregate(t_VAF.75 ~ mc3_exome_barcode, FUN=quantile, 0.75, na.rm=T, data=maf)),
                                                              aggregate(t_VAF.1 ~ mc3_exome_barcode, FUN=quantile, 0.1, na.rm=T, data=maf)),
                                                        aggregate(t_VAF.9 ~ mc3_exome_barcode, FUN=quantile, 0.9, na.rm=T, data=maf)),
                                                  aggregate(p_VAF.25 ~ mc3_exome_barcode, FUN=quantile, 0.25, na.rm=T, data=maf)),
                                            aggregate(p_VAF.75 ~ mc3_exome_barcode, FUN=quantile, 0.75, na.rm=T, data=maf)),
                                      aggregate(p_VAF.1 ~ mc3_exome_barcode, FUN=quantile, 0.1, na.rm=T, data=maf)),
                                aggregate(p_VAF.9 ~ mc3_exome_barcode, FUN=quantile, 0.9, na.rm=T, data=maf)),
                          aggregate(p_VAF ~ mc3_exome_barcode, FUN=mean, data=maf)),
                    aggregate( both ~ mc3_exome_barcode, FUN=mean,  data=maf)),
              aggregate( Start_Position ~ mc3_exome_barcode, FUN=length,  data=maf))

library(ranger)
train_ix = sample(1:nrow(msamp), ceiling(0.6*nrow(msamp)))
strain = msamp[train_ix,]
stest = msamp[-train_ix,]
rf.0 = ranger(both ~ . - mc3_exome_barcode - Start_Position, num.trees =1000, data=strain, importance = "permutation")
p0 = predict(rf.0, data=stest)
stest$predicted = p0$predictions
# plot(stest$predicted, stest$both)

# Predict the fraction of co-called variants by sample from the PCAWG and MC3 VAF distribution in the test set
print(cor(stest$predicted, stest$both, method = "spearman"))

# Out of bag R-squared 
rf.0$r.squared

### 5. GC content
maf$uvar = paste(maf$my.chr, maf$my.pos, maf$my.alt, sep='_')
cadd$uvar = paste(cadd$Chrom, cadd$Pos, cadd$Alt, sep="_")

cadd = cadd[cadd$uvar %in% maf$uvar,]
cadd = cadd[sample(1:nrow(cadd)),]
cadd = cadd[!duplicated(cadd$uvar),]
names(cadd)[c(11,95)] = c("cadd.Consequence", "cadd.CCDS")
maf$unshared = !maf$both
m = merge(cadd,maf)
m$GC.bin = ceiling(m$GC*50)/50
m$yesTcga = m$both | m$t_only
m$yetPcawg = m$both | m$p_only
m$either = T

t.depth.by.GC = aggregate(t_tot ~ GC.bin, FUN=quantile, 0.5, data=m, na.rm=T)
p.depth.by.GC = aggregate(p_tot ~ GC.bin, FUN=quantile, 0.5, data=m, na.rm=T)
n.depth.by.GC = aggregate(n_depth ~ GC.bin, FUN=quantile, 0.5, data=m, na.rm=T)
t.capt.by.GC = aggregate(yesTcga ~ GC.bin, FUN=mean, data=m, na.rm=T)
p.capt.by.GC = aggregate(yetPcawg ~ GC.bin, FUN=mean, data=m, na.rm=T)
p.var.GC = aggregate(either ~ GC.bin, FUN=sum, data=m, na.rm=T)



m.depth.by.GC = merge(merge(merge(merge(merge(t.depth.by.GC, p.depth.by.GC, all=T), n.depth.by.GC, all=T),
                                  t.capt.by.GC, all=T), p.capt.by.GC, all=T), p.var.GC, all=T)


pdf("depthByGCbin.pdf", useDingbats = F)
plot(m.depth.by.GC$GC.bin, m.depth.by.GC$p_tot, ylim=c(0,110), col="red"
     , xlab="GC content", ylab="Median coverage at observed variants")
points(m.depth.by.GC$GC.bin, m.depth.by.GC$t_tot, col="blue")
dev.off()

pdf("RelativeCallRateByGC.pdf", useDingbats = F)
plot(m.depth.by.GC$GC.bin, m.depth.by.GC$yetPcawg, ylim=c(0,1), col="red",
     xlab = "GC content", ylab="Fraction called of merge set")
points(m.depth.by.GC$GC.bin, m.depth.by.GC$yesTcga, col="blue")
dev.off()

# Number of PCAWG variant calls in GC extreme regions
print(sum(m.depth.by.GC$numPcawg[m.depth.by.GC$GC.bin <= 0.3 | m.depth.by.GC$GC.bin >= 0.7]))
# Number of MC3 variant calls in GC extreme regions
print(sum(m.depth.by.GC$numTcga[m.depth.by.GC$GC.bin <= 0.3 | m.depth.by.GC$GC.bin >= 0.7]))
