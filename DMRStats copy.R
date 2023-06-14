setwd("~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/")

CG_1217 <- read.csv("CG_DMRs_2012_2017_0.0001.csv", header = TRUE)
CG_0817 <- read.csv("CG_DMRs_2008_2017_0.0001.csv", header = TRUE)
CG_0812 <- read.csv("CG_DMRs_2008_2012_0.0001.csv", header = TRUE)

CHG_1217 <- read.csv("CHG_DMRs_2012_2017_0.0001.csv", header = TRUE)
CHG_0817 <- read.csv("CHG_DMRs_2008_2017_0.0001.csv", header = TRUE)
CHG_0812 <- read.csv("CHG_DMRs_2008_2012_0.0001.csv", header = TRUE)

CHH_1217 <- read.csv("CHH_2012_2017_0.0001.csv", header = TRUE)
CHH_0817 <- read.csv("CHH_2008_2017_0.0001.csv", header = TRUE)
CHH_0812 <- read.csv("CHH_2008_2012_0.0001.csv", header = TRUE)

summary(CG_0812)
summary(CG_0817)
summary(CG_1217)
str(CG_0812)
str(CG_0817)
str(CG_1217)

nrow(CG_0812[CG_0812$diff.Methy<0,])
nrow(CG_0817[CG_0817$diff.Methy<0,])
nrow(CG_1217[CG_1217$diff.Methy<0,])
nrow(CG_0812[CG_0812$diff.Methy>0,])
nrow(CG_0817[CG_0817$diff.Methy>0,])
nrow(CG_1217[CG_1217$diff.Methy>0,])

nrow(CG_0812[CG_0812$chr == "chr1",])
nrow(CG_0812[CG_0812$chr == "chr2",])
nrow(CG_0812[CG_0812$chr == "chr3",])
nrow(CG_0812[CG_0812$chr == "chr4",])
nrow(CG_0812[CG_0812$chr == "chr5",])
nrow(CG_0812[CG_0812$chr == "chr6",])
nrow(CG_0812[CG_0812$chr == "chr7",])
nrow(CG_0812[CG_0812$chr == "chr8",])

nrow(CG_0817[CG_0817$chr == "chr1",])
nrow(CG_0817[CG_0817$chr == "chr2",])
nrow(CG_0817[CG_0817$chr == "chr3",])
nrow(CG_0817[CG_0817$chr == "chr4",])
nrow(CG_0817[CG_0817$chr == "chr5",])
nrow(CG_0817[CG_0817$chr == "chr6",])
nrow(CG_0817[CG_0817$chr == "chr7",])
nrow(CG_0817[CG_0817$chr == "chr8",])

nrow(CG_1217[CG_1217$chr == "chr1",])
nrow(CG_1217[CG_1217$chr == "chr2",])
nrow(CG_1217[CG_1217$chr == "chr3",])
nrow(CG_1217[CG_1217$chr == "chr4",])
nrow(CG_1217[CG_1217$chr == "chr5",])
nrow(CG_1217[CG_1217$chr == "chr6",])
nrow(CG_1217[CG_1217$chr == "chr7",])
nrow(CG_1217[CG_1217$chr == "chr8",])

mean(CG_0812$length)
mean(CG_0817$length)
mean(CG_1217$length)

summary(CHG_0812)
summary(CHG_0817)
summary(CHG_1217)
str(CHG_0812)
str(CHG_0817)
str(CHG_1217)

nrow(CHG_0812[CHG_0812$diff.Methy<0,])
nrow(CHG_0817[CHG_0817$diff.Methy<0,])
nrow(CHG_1217[CHG_1217$diff.Methy<0,])
nrow(CHG_0812[CHG_0812$diff.Methy>0,])
nrow(CHG_0817[CHG_0817$diff.Methy>0,])
nrow(CHG_1217[CHG_1217$diff.Methy>0,])

nrow(CHG_0812[CHG_0812$chr == "chr1",])
nrow(CHG_0812[CHG_0812$chr == "chr2",])
nrow(CHG_0812[CHG_0812$chr == "chr3",])
nrow(CHG_0812[CHG_0812$chr == "chr4",])
nrow(CHG_0812[CHG_0812$chr == "chr5",])
nrow(CHG_0812[CHG_0812$chr == "chr6",])
nrow(CHG_0812[CHG_0812$chr == "chr7",])
nrow(CHG_0812[CHG_0812$chr == "chr8",])

nrow(CHG_0817[CHG_0817$chr == "chr1",])
nrow(CHG_0817[CHG_0817$chr == "chr2",])
nrow(CHG_0817[CHG_0817$chr == "chr3",])
nrow(CHG_0817[CHG_0817$chr == "chr4",])
nrow(CHG_0817[CHG_0817$chr == "chr5",])
nrow(CHG_0817[CHG_0817$chr == "chr6",])
nrow(CHG_0817[CHG_0817$chr == "chr7",])
nrow(CHG_0817[CHG_0817$chr == "chr8",])

nrow(CHG_1217[CHG_1217$chr == "chr1",])
nrow(CHG_1217[CHG_1217$chr == "chr2",])
nrow(CHG_1217[CHG_1217$chr == "chr3",])
nrow(CHG_1217[CHG_1217$chr == "chr4",])
nrow(CHG_1217[CHG_1217$chr == "chr5",])
nrow(CHG_1217[CHG_1217$chr == "chr6",])
nrow(CHG_1217[CHG_1217$chr == "chr7",])
nrow(CHG_1217[CHG_1217$chr == "chr8",])


mean(CHG_0812$length)
mean(CHG_0817$length)
mean(CHG_1217$length)

summary(CHH_0812)
summary(CHH_0817)
summary(CHH_1217)
str(CHH_0812)
str(CHH_0817)
str(CHH_1217)

nrow(CHH_0812[CHH_0812$diff.Methy<0,])
nrow(CHH_0817[CHH_0817$diff.Methy<0,])
nrow(CHH_1217[CHH_1217$diff.Methy<0,])
nrow(CHH_0812[CHH_0812$diff.Methy>0,])
nrow(CHH_0817[CHH_0817$diff.Methy>0,])
nrow(CHH_1217[CHH_1217$diff.Methy>0,])

nrow(CHH_0812[CHH_0812$chr == "chr1",])
nrow(CHH_0812[CHH_0812$chr == "chr2",])
nrow(CHH_0812[CHH_0812$chr == "chr3",])
nrow(CHH_0812[CHH_0812$chr == "chr4",])
nrow(CHH_0812[CHH_0812$chr == "chr5",])
nrow(CHH_0812[CHH_0812$chr == "chr6",])
nrow(CHH_0812[CHH_0812$chr == "chr7",])
nrow(CHH_0812[CHH_0812$chr == "chr8",])

nrow(CHH_0817[CHH_0817$chr == "chr1",])
nrow(CHH_0817[CHH_0817$chr == "chr2",])
nrow(CHH_0817[CHH_0817$chr == "chr3",])
nrow(CHH_0817[CHH_0817$chr == "chr4",])
nrow(CHH_0817[CHH_0817$chr == "chr5",])
nrow(CHH_0817[CHH_0817$chr == "chr6",])
nrow(CHH_0817[CHH_0817$chr == "chr7",])
nrow(CHH_0817[CHH_0817$chr == "chr8",])

nrow(CHH_1217[CHH_1217$chr == "chr1",])
nrow(CHH_1217[CHH_1217$chr == "chr2",])
nrow(CHH_1217[CHH_1217$chr == "chr3",])
nrow(CHH_1217[CHH_1217$chr == "chr4",])
nrow(CHH_1217[CHH_1217$chr == "chr5",])
nrow(CHH_1217[CHH_1217$chr == "chr6",])
nrow(CHH_1217[CHH_1217$chr == "chr7",])
nrow(CHH_1217[CHH_1217$chr == "chr8",])


mean(CHH_0812$length)
mean(CHH_0817$length)
mean(CHH_1217$length)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

getmode(CG_0812$length)
getmode(CG_0817$length)
getmode(CG_1217$length)

getmode(CHG_0812$length)
getmode(CHG_0817$length)
getmode(CHG_1217$length)

getmode(CHH_0812$length)
getmode(CHH_0817$length)
getmode(CHH_1217$length)

CG_DMRdist <- read.csv("CG_DMRDist_Chrm.csv", header = TRUE)
str(CG_DMRdist)
CG_DMRdist$Chromsome <- as.factor(CG_DMRdist$Chromsome)
CG_DMRdist$Comparison <- as.factor(CG_DMRdist$Comparison)

library(ggplot2)
cbp1 <- c("#D55E00", "#999999", "#E69F00", "#0072B2", "#CC79A7", "#56B4E9", "#009E73","#F0E442")

CG_Dist <- ggplot(CG_DMRdist, aes(x = Chromsome, y = Num_DMRs)) + geom_point(aes(color = Comparison, shape = Comparison), size=3) + labs(x = "", y = "Number of DMRs") + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + scale_color_manual(values = cbp1)
pdf("DMR_CGDist_Chrom.pdf")
CG_Dist
dev.off()
CHG_DMRdist <- read.csv("CHG_DMRDist_Chrm.csv", header = TRUE)
str(CHG_DMRdist)
CHG_DMRdist$Chromsome <- as.factor(CHG_DMRdist$Chromsome)
CHG_DMRdist$Comparison <- as.factor(CHG_DMRdist$Comparison)

CHG_Dist <- ggplot(CHG_DMRdist, aes(x = Chromsome, y = Num_DMRs)) + geom_point(aes(color = Comparison, shape = Comparison), size=3) + labs(x = "", y = "Number of DMRs") + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + scale_color_manual(values = cbp1)

CHH_DMRdist <- read.csv("CHH_DMRDist_Chrm.csv", header = TRUE)
str(CHH_DMRdist)
CHH_DMRdist$Chromsome <- as.factor(CHH_DMRdist$Chromsome)
CHH_DMRdist$Comparison <- as.factor(CHH_DMRdist$Comparison)

CHH_Dist <- ggplot(CHH_DMRdist, aes(x = Chromsome, y = Num_DMRs)) + geom_point(aes(color = Comparison, shape = Comparison), size=3) + labs(x = "", y = "Number of DMRs") + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")+ scale_color_manual(values = cbp1)

library(multipanelfigure)
figureDist <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figureDist %<>% fill_panel(CG_Dist)
figureDist %<>% fill_panel(CHG_Dist)
figureDist %<>% fill_panel(CHH_Dist)
figureDist
pdf("DMRDist_Chromosome.pdf")
figureDist
dev.off()
