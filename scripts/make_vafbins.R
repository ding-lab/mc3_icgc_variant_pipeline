library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)


args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Full.tsv and output should be supplied (input file).n", call.=FALSE)
}

data = fread(input=args[1], sep="\t", header=TRUE, na.strings="NA",colClasses=list(character=c(128,135,139)))
data$match = ifelse(!is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$PCAWG_only = ifelse(is.na(data$Chromosome) & !is.na(data$"Chromosome:1"),1,0)
data$MC3_only = ifelse(is.na(data$"Chromosome:1") & !is.na(data$Chromosome),1,0)
data$i_VAFfixed = sapply(strsplit(as.character(data$i_VAF), "|", fixed=T), function(x) mean(as.numeric(x)))



##### TRY TO GET MORE GRANULARITY OUT OF THIS FIGURE 
PASS <- data
PASS$binned_iVAF = ifelse(is.na(PASS$i_VAFfixed),"NA",ceiling(as.numeric(PASS$i_VAFfixed)*20) / 20)


PASS$m_VAF = ifelse(is.na(PASS$t_depth),"NA",PASS$t_alt_count/PASS$t_depth)
PASS$binned_mVAF = ifelse(is.na(PASS$m_VAF),"NA",ceiling(as.numeric(PASS$m_VAF)*20) / 20)

PASS$agreement = ifelse(PASS$binned_iVAF == PASS$binned_mVAF, 1, 0)

combinedVAF <- function(vec){
    a = sum(na.omit(as.numeric(c(vec['t_alt_count'],vec['t_alt_count:1']))))
    b = sum(na.omit(as.numeric(c(vec['t_depth'],vec['t_alt_count:1'],vec['t_ref_count:1']))))
    return (a/b)
}

PASS$cVAF = apply(PASS,1,combinedVAF)
PASS$binned_cVAF = ifelse(is.na(PASS$cVAF),"NA",ceiling(as.numeric(PASS$cVAF)*20) / 20)


PASS$Matched = ifelse(PASS$match == 1, "Match", ifelse(PASS$MC3_only, "MC3_uniq", "PCAWG_uniq"))

dat3 <- PASS %>%
  group_by(binned_cVAF) %>% mutate(n = n()) %>%
  mutate(label = paste0(binned_cVAF,' N=',n))

dat3$n <- NULL


summary2 = dat3 %>% group_by(Matched,label) %>%
  tally %>%
  group_by(label) %>%
  mutate(pct = n/sum(n),
         n.pos = cumsum(n) - 0.5)


#p <- ggplot(summary[which(summary$label != "NA\nN = 59232"),], aes(x=label, y=pct, fill=TCGA_Uniq))
p <- ggplot(summary2, aes(x=label, y=pct, fill=Matched))
p <- p+geom_bar(stat="identity")
#p <- p+geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%"), y=ifelse(match,pct*.5,.94),angle=90),colour="white")
p <- p+labs(y="Percent of mutations",x="Combined MC3+PCAWG VAF binned by 5% intervals")
p <- p+scale_y_continuous(labels=percent)
p <- p+scale_fill_manual(values=c("#4DAF4A","#5981A6","#BF312E"))
p <- p+theme_minimal()
p <- p+theme(legend.position="bottom",axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
p



pdf(args[2],height=4,width=4,useDingbats=F)
print(p)
dev.off()

