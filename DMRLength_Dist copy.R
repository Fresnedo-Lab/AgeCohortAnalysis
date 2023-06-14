setwd("~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/")

CG_0812 <- read.csv("CG_DMRs_2008_2012_0.0001.csv", header = TRUE)
CG_0817 <- read.csv("CG_DMRs_2008_2017_0.0001.csv", header = TRUE)
CG_1217 <- read.csv("CG_DMRs_2012_2017_0.0001.csv", header = TRUE)

CHG_0812 <- read.csv("CHG_DMRs_2008_2012_0.0001.csv", header = TRUE)
CHG_0817 <- read.csv("CHG_DMRs_2008_2017_0.0001.csv", header = TRUE)
CHG_1217 <- read.csv("CHG_DMRs_2012_2017_0.0001.csv", header = TRUE)

CHH_0812 <- read.csv("CHH_2008_2012_0.0001.csv", header = TRUE)
CHH_0817 <- read.csv("CHH_2008_2017_0.0001.csv", header = TRUE)
CHH_1217 <- read.csv("CHH_2012_2017_0.0001.csv", header = TRUE)

max(CG_0812$length)
max(CG_0817$length)
max(CG_1217$length)

max(CHG_0812$length)
max(CHG_0817$length)
max(CHG_1217$length)

max(CHH_0812$length)
max(CHH_0817$length)
max(CHH_1217$length)

min(CG_0812$length)
min(CG_0817$length)
min(CG_1217$length)

min(CHG_0812$length)
min(CHG_0817$length)
min(CHG_1217$length)

min(CHH_0812$length)
min(CHH_0817$length)
min(CHH_1217$length)

mean(CG_0812$length)
mean(CG_0817$length)
mean(CG_1217$length)

mean(CHG_0812$length)
mean(CHG_0817$length)
mean(CHG_1217$length)

mean(CHH_0812$length)
mean(CHH_0817$length)
mean(CHH_1217$length)

median(CG_0812$length)
median(CG_0817$length)
median(CG_1217$length)

median(CHG_0812$length)
median(CHG_0817$length)
median(CHG_1217$length)

median(CHH_0812$length)
median(CHH_0817$length)
median(CHH_1217$length)

getmode(CG_0812$length)
getmode(CG_0817$length)
getmode(CG_1217$length)

getmode(CHG_0812$length)
getmode(CHG_0817$length)
getmode(CHG_1217$length)

getmode(CHH_0812$length)
getmode(CHH_0817$length)
getmode(CHH_1217$length)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


library(ggplot2)
library(patchwork)

CG_0812_P <- ggplot(CG_0812, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CG 11-7",x="", y = "")
CG_0817_P <- ggplot(CG_0817, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CG 11-2",x="", y = "Count")
CG_1217_P <- ggplot(CG_1217, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CG 7-2",x="", y = "")

CHG_0812_P <- ggplot(CHG_0812, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHG 11-7",x="", y = "")
CHG_0817_P <- ggplot(CHG_0817, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHG 11-2",x="", y = "Count")
CHG_1217_P <- ggplot(CHG_1217, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHG 7-2",x="", y = "")

CHH_0812_P <- ggplot(CHH_0812, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHH 11-7",x="DMR Length (base pairs)", y = "")
CHH_0817_P <- ggplot(CHH_0817, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHH 11-2",x="", y = "Count")
CHH_1217_P <- ggplot(CHH_1217, aes(x=length)) + geom_histogram(binwidth = 15, colour="black", fill="white") + theme_bw() + labs(title="CHH 7-2",x="", y = "")

patchwork_CG <- CG_0817_P + CG_0812_P + CG_1217_P
pdf("CG_DMRDistribution.pdf")
patchwork_CG + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()

patchwork_CHG <- CHG_0817_P + CHG_0812_P + CHG_1217_P
pdf("CHG_DMRDistribution.pdf")
patchwork_CHG + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()

patchwork_CHH <- CHH_0817_P + CHH_0812_P + CHH_1217_P
pdf("CHH_DMRDistribution.pdf")
patchwork_CHH + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()

patchwork_All <- CG_0817_P + CG_0812_P + CG_1217_P + CHG_0817_P + CHG_0812_P + CHG_1217_P + CHH_0817_P + CHH_0812_P + CHH_1217_P
pdf("DMRDistribution_All.pdf")
patchwork_All+plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12, face = "bold"))
dev.off()
