library(data.table)
library(ggplot2)
library(dplyr)


tcga_colors <- c(
  `BRCA`  = "#ED2891",
  `PRAD`  = "#7E1918",
  `TGCT`  = "#BE1E2D",
  `KICH`  = "#ED1C24",
  `KIRP`  = "#EA7075",
  `KIRC`  = "#F8AFB3",
  `BLCA`  = "#FAD2D9",
  `OV`    = "#D97D25",
  `UCS`   = "#F89420",
  `CESC`  = "#F6B667",
  `UCEC`  = "#FBE3C7",
  `THCA`  = "#F9ED32",
  `PCPG`  = "#E8C51D",
  `ACC`   = "#C1A72F",
  `SKCM`  = "#BBD642",
  `UVM`   = "#009444",
  `HNSC`  = "#97D1A9",
  `SARC`  = "#00A99D",
  `ESCA`  = "#007EB5",
  `STAD`  = "#00AEEF",
  `COAD`  = "#9EDDF9",
  `READ`  = "#DAF1FC",
  `CHOL`  = "#104A7F",
  `PAAD`  = "#6E7BA2",
  `LIHC`  = "#CACCDB",
  `MESO`  = "#542C88",
  `LUSC`  = "#A084BD",
  `LUAD`  = "#D3C3E0",
  `GBM`   = "#B2509E",
  `LGG`   = "#D49DC7",
  `DLBC`  = "#3953A4",
  `LAML`  = "#754C29",
  `THYM`  = "#CEAC8F")



tcga_cols <- function(...) {
  cols <- c(...)



  tcga_colors[cols]
}

tcga_palettes <- list(
  `main`  = tcga_cols("KICH", "STAD", "MESO"),

  `cool`  = tcga_cols("SARC", "ESCA", "STAD", "COAD", "READ","CHOL","PAAD","LIHC","MESO","LUSC"),

  `hot`   = tcga_cols("BRCA", "PRAD", "TCGT","KICH","KIRP","OV","UCS","CESC","THCA","PCPG"),

  `mixed` = tcga_cols("BRCA", "PRAD", "THCA", "SKCM", "UVM","SARC","ESCA","CHOL","MESO","GBM","DLBC"),

  `neutral`  = tcga_cols("LIHC", "ACC","LAML","THYM","UCEC")
)

tcga_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- tcga_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  colorRampPalette(pal, ...)
}


#Accessing these palettes
tcga_pal("cool")
tcga_pal("cool")(10)

#SET THIS UP for GGPLOT2 -- TCGA
scale_color_tcga <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- tcga_pal(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("colour", paste0("tcga_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}



scale_fill_tcga <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- tcga_pal(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("fill", paste0("tcga_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}
#For references this code was taken and adapted from here:
#https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2



#columns that I want to keep and compare
#histological_type
#race
#tumor_tissue_site
#bcr_patient_barcode
#acronym
#gender

#Functions 
'%!in%' <- function(x,y)!('%in%'(x,y))


#Arugments
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("Something went wrong with the snakemake call and there were not enough files provided", call.=FALSE)
}


clin <- fread(args[1],header=T)
my_clin <- data.frame("ID"=clin$bcr_patient_barcode,"RACE"=clin$race,"TT"=clin$tumor_tissue_site,"CODE"=clin$acronym,"GENDER"=clin$gender,"HIST"=clin$histological_type)

WGS <- fread(args[2],header=T)
wgschar12 <- substr(WGS$mc3_exome_barcode,1,12)


my_clin$InWGS <- ifelse(my_clin$ID %in% wgschar12, "YES", "NO")


cantypes <- unique(my_clin$CODE)


#For testing 
#can = "BRCA"
#hist = "Infiltrating Ductal Carcinoma"

#For each cancer 
diffs <- NULL
for(can in cantypes){
    TEST <- my_clin[which(my_clin$CODE == can),]
    can_tot <- dim(TEST)[1]
    wgs_tot <- dim(TEST[which(TEST$InWGS == "YES"),])[1]
    hists <- unique(TEST$HIST)
    #For each histological type
    for(hist in hists){
        TESTa <- TEST[which(TEST$HIST == hist),]
        hc_tot <- dim(TESTa)[1]
        hw_tot <- dim(TESTa[which(TESTa$InWGS == "YES"),])[1]
        a = c(hw_tot, (wgs_tot-hw_tot))
        b = c(hc_tot, (can_tot-hc_tot))
        d = rbind(a,b)
        myf <- fisher.test(d)
        o = data.frame("Cancer"=can,"HISTOLOGY"=hist,"Pval"=myf$p.value,"OddsRat"=myf$estimate,"Conf95Low"=myf$conf.int[1],"Conf95High"=myf$conf.int[2],"WGSratio"=a[1]/sum(a),"TCGAratio"=b[1]/sum(b),"WGS_yes"=a[1],"WGS_tot"=wgs_tot,"TCGA_yes"=b[1],"TCGA_tot"=can_tot)
        diffs <- rbind(diffs,o)
    }
}




p <- ggplot(diffs,aes(x=log2(OddsRat),y=-log10(Pval),color=Cancer))
p <- p + geom_point(alpha=.85,size=4.5,shape=16)
p <- p + geom_hline(yintercept=-log10(0.05),color="red")
p <- p + geom_vline(xintercept=0,color="red")
#p <- p + discrete_color_tcga(tcga_cls)
p <- p + scale_color_manual(values=tcga_colors)
p <- p + geom_text(aes(label=ifelse(-log10(Pval)> 1.3, as.character(paste(Cancer,HISTOLOGY,sep=": ")),'')),color="black",hjust=0,vjust=.5,angle = 90)
p <- p + theme_minimal()
p <- p + theme(legend.position = "none", legend.box = "horizontal")
p <- p + xlab("log2(Odds Ratio)") + ylab("-log10(P-value)")
p

pdf(args[3],height=3,width=5,useDingbats=F)
print(p)
dev.off()




p <- ggplot(diffs,aes(x=log2(OddsRat),y=-log10(Pval),color=Cancer))
p <- p + geom_point(alpha=.5,size=4.5)
p <- p + geom_hline(yintercept=-log10(0.05),color="red")
#p <- p + discrete_color_tcga(tcga_cls)
p <- p + scale_color_manual(values=tcga_colors)
p <- p + geom_text(aes(label=ifelse(-log10(Pval)> 1.3, as.character(HISTOLOGY),'')),color="black",hjust=0,vjust=0,angle = 90)
p <- p + theme_minimal()
p <- p + theme(legend.position = "bottom", legend.box = "horizontal")
p

pdf(args[4]"SelectionBias.v1.withlegend.pdf",height=4,width=6,useDingbats=F)
print(p)
dev.off()




#Here is a histogram of all of the TCGA cancter types and where they fit into the scheme of things 

t <- ggplot(mi_counts,aes(x=Cancer,y=sum.WGS_yes.,fill=Cancer))
t <- t+geom_bar(stat="identity")
t <- t+geom_text(stat='identity', aes(label=sum.WGS_yes.), vjust=-1)
t <- t+theme_classic()
t <- t+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5),legend.position="None")
t <- t+scale_fill_manual(values=tcga_colors)
t <- t+ylab("Number of samples")+xlab("TCGA Cancer types")
t

pdf(args[5],height=3,width=8,useDingbats=F)
print(t)
dev.off()

#Write out the datatable 

write.table(diffs,file=args[6],sep="\t", quote=F, row.names=F)





