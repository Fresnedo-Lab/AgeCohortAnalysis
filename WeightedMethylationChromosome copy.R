library(multipanelfigure)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)
library(wesanderson)


setwd("~/Documents/Almond/AgeCohortProject/")

CG_chr <- read.delim("CG_wide.txt", header = TRUE)
str(CG_chr)
CG_chr$Age <- as.factor(CG_chr$Age)
CG_chr <- CG_chr[-2]
melt(CG_chr, variable.name = "chromosome", value.name = "per_meth")
reCG_chr <- melt(CG_chr, variable.name = "chromosome", value.name = "per_meth")

cbp1 <- c("#D55E00", "#999999", "#E69F00", "#0072B2", "#CC79A7", "#56B4E9", "#009E73","#F0E442")

pdf("WeightedMeth_Chromosome_CG.pdf")
CGWmeth <- ggplot(reCG_chr, aes(x = chromosome, y = per_meth, fill = Age)) + geom_boxplot() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + scale_fill_manual(values = cbp1) + ylab("Proportation Weighted CG Methylation")
dev.off()
gy_logit_CG <- betareg(per_meth~Age*chromosome, data = reCG_chr)
ggplot(reCG_chr, aes(x=chromosome, y=per_meth)) + geom_boxplot()
summary(gy_logit_CG)
marginal_CG <- emmeans(gy_logit_CG, ~Age*chromosome)
pairs(marginal_CG)
pwpp(marginal_CG, by = "chromosome")

CHG_chr <- read.delim("CHG_wide.txt", header = TRUE)
str(CHG_chr)
CHG_chr$Age <- as.factor(CHG_chr$Age)
CHG_chr <- CHG_chr[-2]
melt(CHG_chr, variable.name = "chromosome", value.name = "per_meth")
reCHG_chr <- melt(CHG_chr, variable.name = "chromosome", value.name = "per_meth")
pdf("WeightedMeth_Chromosome_CHG.pdf")
CHGWmeth<- ggplot(reCHG_chr, aes(x = chromosome, y = per_meth, fill = Age)) + geom_boxplot() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + scale_fill_manual(values = cbp1) + ylab("Proportation Weighted CHG Methylation")
dev.off()
gy_logit_CHG <- betareg(per_meth~Age*chromosome, data = reCHG_chr)
ggplot(reCHG_chr, aes(x=chromosome, y=per_meth)) + geom_boxplot()

summary(gy_logit_CHG)
marginal_CHG <- emmeans(gy_logit_CHG, ~Age*chromosome)
pairs(marginal_CHG)
pwpp(marginal_CHG, by = "chromosome")
pwpp(marginal_CHG)


CHH_chr <- read.delim("CHH_wide.txt", header = TRUE)
str(CHH_chr)
CHH_chr$Age <- as.factor(CHH_chr$Age)
CHH_chr <- CHH_chr[-2]
melt(CHH_chr, variable.name = "chromosome", value.name = "per_meth")
reCHH_chr <- melt(CHH_chr, variable.name = "chromosome", value.name = "per_meth")
redCHH_chr <- CHH_chr[c(-21,-30),]
redCHH_chr <- melt(redCHH_chr, variable.name = "chromosome", value.name = "per_meth")
pdf("WeightedMeth_Chromosome_CHH.pdf")
CHHWmeth<- ggplot(reCHH_chr, aes(x = chromosome, y = per_meth, fill = Age)) + geom_boxplot() + theme_bw() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + scale_fill_manual(values = cbp1) + ylab("Proportation Weighted CHH Methylation")
dev.off()
gy_logit_CHH <- betareg(per_meth~Age*chromosome, data = reCHH_chr)
ggplot(reCHH_chr, aes(x=chromosome, y=per_meth)) + geom_boxplot()

ggplot(redCHH_chr, aes(x = chromosome, y = per_meth, fill = Age)) + geom_boxplot() + theme(axis.title.x = element_blank()) + scale_fill_manual(values=wes_palette(n=3, name="Royal2")) + ylab("Proportation Weighted CHH Methylation")
ggplot(redCHH_chr, aes(x=chromosome, y=per_meth)) + geom_boxplot()

summary(gy_logit_CHH)
marginal_CHH <- emmeans(gy_logit_CHH, ~Age*chromosome)
pairs(marginal_CHH)
pwpp(marginal_CHH, by = "chromosome")

figureChr <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figureChr %<>% fill_panel(CGWmeth)
figureChr %<>% fill_panel(CHGWmeth)
figureChr %<>% fill_panel(CHHWmeth)
figureChr
pdf("PercentMethylation_ChromosomePlot.pdf")
figureChr
dev.off()



stat.test1 <- tibble::tribble(
  ~group1, ~group2,   ~p.value,
  "2",     "7", 0.99,
  "2",     "11", 0.01,
  "7",     "11", 0.01
)
stat.test1
plot1 <- ggplot(CG_chr, aes(x = Age, y = Chr1)) + geom_boxplot() + stat_pvalue_manual(stat.test, label = "p.value", y.position = 0.6, step.increase = 0.1)
plot(pairs(marginal1)) + theme_bw() + 
  labs(y="", x = "Estimated mean difference")
plot(marginal1, comparisons = TRUE) + theme_bw() + 
  labs(y = "", x = "Estimated marginal mean (Weighted CG methylation Chr1)")
plot2 <- ggplot(CG_chr, aes(x = Age, y = Chr2)) + geom_boxplot() 
plot3 <- ggplot(CG_chr, aes(x = Age, y = Chr3)) + geom_boxplot() 
plot4 <- ggplot(CG_chr, aes(x = Age, y = Chr4)) + geom_boxplot() 
plot5 <- ggplot(CG_chr, aes(x = Age, y = Chr5)) + geom_boxplot() 
plot6 <- ggplot(CG_chr, aes(x = Age, y = Chr6)) + geom_boxplot() 
plot7 <- ggplot(CG_chr, aes(x = Age, y = Chr7)) + geom_boxplot() 
plot8 <- ggplot(CG_chr, aes(x = Age, y = Chr8)) + geom_boxplot() 


figure1a <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure1a %<>% fill_panel(plot1)
figure1a %<>% fill_panel(plot2)
figure1a %<>% fill_panel(plot3)
figure1a %<>% fill_panel(plot4)
figure1a
pdf("PercentMethylation_Chr1_4_CG.pdf")
figure1a
dev.off()

figure1b <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure1b %<>% fill_panel(plot5)
figure1b %<>% fill_panel(plot6)
figure1b %<>% fill_panel(plot7)
figure1b %<>% fill_panel(plot8)
figure1b
pdf("PercentMethylation_Chr5_8_CG.pdf")
figure1b
dev.off()

plot9 <- ggplot(CHG_chr, aes(x = Age, y = Chr1)) + geom_boxplot() 
plot10 <- ggplot(CHG_chr, aes(x = Age, y = Chr2)) + geom_boxplot() 
plot11 <- ggplot(CHG_chr, aes(x = Age, y = Chr3)) + geom_boxplot() 
plot12 <- ggplot(CHG_chr, aes(x = Age, y = Chr4)) + geom_boxplot() 
plot13 <- ggplot(CHG_chr, aes(x = Age, y = Chr5)) + geom_boxplot() 
plot14 <- ggplot(CHG_chr, aes(x = Age, y = Chr6)) + geom_boxplot() 
plot15 <- ggplot(CHG_chr, aes(x = Age, y = Chr7)) + geom_boxplot() 
plot16 <- ggplot(CHG_chr, aes(x = Age, y = Chr8)) + geom_boxplot() 


figure2a <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure2a %<>% fill_panel(plot9)
figure2a %<>% fill_panel(plot10)
figure2a %<>% fill_panel(plot11)
figure2a %<>% fill_panel(plot12)
figure2a
pdf("PercentMethylation_Chr1_4_CHG.pdf")
figure2a
dev.off()

figure2b <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure2b %<>% fill_panel(plot13)
figure2b %<>% fill_panel(plot14)
figure2b %<>% fill_panel(plot15)
figure2b %<>% fill_panel(plot16)
figure2b
pdf("PercentMethylation_Chr5_8_CHG.pdf")
figure2b
dev.off()

plot17 <- ggplot(CHH_chr, aes(x = Age, y = Chr1)) + geom_boxplot() 
plot18 <- ggplot(CHH_chr, aes(x = Age, y = Chr2)) + geom_boxplot() 
plot19 <- ggplot(CHH_chr, aes(x = Age, y = Chr3)) + geom_boxplot() 
plot20 <- ggplot(CHH_chr, aes(x = Age, y = Chr4)) + geom_boxplot() 
plot21 <- ggplot(CHH_chr, aes(x = Age, y = Chr5)) + geom_boxplot() 
plot22 <- ggplot(CHH_chr, aes(x = Age, y = Chr6)) + geom_boxplot() 
plot23 <- ggplot(CHH_chr, aes(x = Age, y = Chr7)) + geom_boxplot() 
plot24 <- ggplot(CHH_chr, aes(x = Age, y = Chr8)) + geom_boxplot() 


figure3a <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure3a %<>% fill_panel(plot17)
figure3a %<>% fill_panel(plot18)
figure3a %<>% fill_panel(plot19)
figure3a %<>% fill_panel(plot20)
figure3a
pdf("PercentMethylation_Chr1_4_CHH.pdf")
figure3a
dev.off()

figure3b <- multi_panel_figure(width = 180, height = 180, columns = 2, rows = 2)
figure3b %<>% fill_panel(plot21)
figure3b %<>% fill_panel(plot22)
figure3b %<>% fill_panel(plot23)
figure3b %<>% fill_panel(plot24)
figure3b
pdf("PercentMethylation_Chr5_8_CHH.pdf")
figure3b
dev.off()

library(betareg)
library(emmeans)
library(lmtest)
library(multcompView)
library(multcomp)

gy_logit1 <- betareg(Chr1~Age, data = CG_chr)
gy_logit2 <- betareg(Chr2~Age, data = CG_chr)
gy_logit3 <- betareg(Chr3~Age, data = CG_chr)
gy_logit4 <- betareg(Chr4~Age, data = CG_chr)
gy_logit5 <- betareg(Chr5~Age, data = CG_chr)
gy_logit6 <- betareg(Chr6~Age, data = CG_chr)
gy_logit7 <- betareg(Chr7~Age, data = CG_chr)
gy_logit8 <- betareg(Chr8~Age, data = CG_chr)

summary(gy_logit1)
summary(gy_logit2)
summary(gy_logit3)
summary(gy_logit4)
summary(gy_logit5)
summary(gy_logit6)
summary(gy_logit7)
summary(gy_logit8)

marginal1 <- emmeans(gy_logit1, ~Age)
marginal2 <- emmeans(gy_logit2, ~Age)
marginal3 <- emmeans(gy_logit3, ~Age)
marginal4 <- emmeans(gy_logit4, ~Age)
marginal5 <- emmeans(gy_logit5, ~Age)
marginal6 <- emmeans(gy_logit6, ~Age)
marginal7 <- emmeans(gy_logit7, ~Age)
marginal8 <- emmeans(gy_logit8, ~Age)

pwpp(marginal1, method = "pairwise")
pwpp(marginal2, method = "pairwise")
pwpp(marginal3, method = "pairwise")
pwpp(marginal4, method = "pairwise")
pwpp(marginal5, method = "pairwise")
pwpp(marginal6, method = "pairwise")
pwpp(marginal7, method = "pairwise")
pwpp(marginal8, method = "pairwise")

pairs(marginal1)
pairs(marginal2)
pairs(marginal3)
pairs(marginal4)
pairs(marginal5)
pairs(marginal6)
pairs(marginal7)
pairs(marginal8)

Sum1 <- cld(marginal1, alpha = 0.05, Letters = letters)
Sum1
Sum2 <- cld(marginal2, alpha = 0.05, Letters = letters)
Sum2
Sum3 <- cld(marginal3, alpha = 0.05, Letters = letters)
Sum3
Sum4 <- cld(marginal4, alpha = 0.05, Letters = letters)
Sum4
Sum5 <- cld(marginal5, alpha = 0.05, Letters = letters)
Sum5
Sum6 <- cld(marginal6, alpha = 0.05, Letters = letters)
Sum6
Sum7 <- cld(marginal7, alpha = 0.05, Letters = letters)
Sum7
Sum8 <- cld(marginal8, alpha = 0.05, Letters = letters)
Sum8

gy_logit9 <- betareg(Chr1~Age, data = CHG_chr)
gy_logit10 <- betareg(Chr2~Age, data = CHG_chr)
gy_logit11 <- betareg(Chr3~Age, data = CHG_chr)
gy_logit12 <- betareg(Chr4~Age, data = CHG_chr)
gy_logit13 <- betareg(Chr5~Age, data = CHG_chr)
gy_logit14 <- betareg(Chr6~Age, data = CHG_chr)
gy_logit15 <- betareg(Chr7~Age, data = CHG_chr)
gy_logit16 <- betareg(Chr8~Age, data = CHG_chr)

summary(gy_logit9)
summary(gy_logit10)
summary(gy_logit11)
summary(gy_logit12)
summary(gy_logit13)
summary(gy_logit14)
summary(gy_logit15)
summary(gy_logit16)

marginal9 <- emmeans(gy_logit9, ~Age)
marginal10 <- emmeans(gy_logit10, ~Age)
marginal11 <- emmeans(gy_logit11, ~Age)
marginal12 <- emmeans(gy_logit12, ~Age)
marginal13 <- emmeans(gy_logit13, ~Age)
marginal14 <- emmeans(gy_logit14, ~Age)
marginal15 <- emmeans(gy_logit15, ~Age)
marginal16 <- emmeans(gy_logit16, ~Age)

pairs(marginal9)
pairs(marginal10)
pairs(marginal11)
pairs(marginal12)
pairs(marginal13)
pairs(marginal14)
pairs(marginal15)
pairs(marginal16)

gy_logit17 <- betareg(Chr1~Age, data = CHH_chr)
gy_logit18 <- betareg(Chr2~Age, data = CHH_chr)
gy_logit19 <- betareg(Chr3~Age, data = CHH_chr)
gy_logit20 <- betareg(Chr4~Age, data = CHH_chr)
gy_logit21 <- betareg(Chr5~Age, data = CHH_chr)
gy_logit22 <- betareg(Chr6~Age, data = CHH_chr)
gy_logit23 <- betareg(Chr7~Age, data = CHH_chr)
gy_logit24 <- betareg(Chr8~Age, data = CHH_chr)

summary(gy_logit17)
summary(gy_logit18)
summary(gy_logit19)
summary(gy_logit20)
summary(gy_logit21)
summary(gy_logit22)
summary(gy_logit23)
summary(gy_logit24)

marginal17 <- emmeans(gy_logit17, ~Age)
marginal18 <- emmeans(gy_logit18, ~Age)
marginal19 <- emmeans(gy_logit19, ~Age)
marginal20 <- emmeans(gy_logit20, ~Age)
marginal21 <- emmeans(gy_logit21, ~Age)
marginal22 <- emmeans(gy_logit22, ~Age)
marginal23 <- emmeans(gy_logit23, ~Age)
marginal24 <- emmeans(gy_logit24, ~Age)

pairs(marginal17)
pairs(marginal18)
pairs(marginal19)
pairs(marginal20)
pairs(marginal21)
pairs(marginal22)
pairs(marginal23)
pairs(marginal24)

Sum1 <- cld(marginal1, alpha = 0.05, Letters = letters)
Sum1
Sum2 <- cld(marginal2, alpha = 0.05, Letters = letters)
Sum2
Sum3 <- cld(marginal3, alpha = 0.05, Letters = letters)
Sum3
Sum4 <- cld(marginal4, alpha = 0.05, Letters = letters)
Sum4
Sum5 <- cld(marginal5, alpha = 0.05, Letters = letters)
Sum5
Sum6 <- cld(marginal6, alpha = 0.05, Letters = letters)
Sum6
Sum7 <- cld(marginal7, alpha = 0.05, Letters = letters)
Sum7
Sum8 <- cld(marginal8, alpha = 0.05, Letters = letters)
Sum8

plot1 <- ggplot(Sum1,                
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr1") +
  
  annotate("text",
           x = 1:3,
           y = Sum1$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum1$.group))

plot2 <- ggplot(Sum2,                ### The data frame to use.
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr2") +
  
  annotate("text",
           x = 1:3,
           y = Sum2$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum2$.group))

plot3 <- ggplot(Sum3,                ### The data frame to use.
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr3") +
  
  annotate("text",
           x = 1:3,
           y = Sum3$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum3$.group))

plot4 <- ggplot(Sum4,                ### The data frame to use.
aes(x = Age,
y = emmean)) +
geom_errorbar(aes(ymin = asymp.LCL,
ymax = asymp.UCL),
width = 0.05,
size  = 0.5) +
geom_point(shape = 15,
size  = 4) +
theme_bw() +
theme(axis.title   = element_text(face  = "bold")) +
ylab("Proportion of CG methylation Chr4") +

annotate("text",
x = 1:3,
 y = Sum4$asymp.UCL + 0.008,
label = gsub(" ", "", Sum4$.group))

plot5 <- ggplot(Sum5,                
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr5") +
  
  annotate("text",
           x = 1:3,
           y = Sum1$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum1$.group))

plot6 <- ggplot(Sum6,                ### The data frame to use.
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr6") +
  
  annotate("text",
           x = 1:3,
           y = Sum2$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum2$.group))

plot7 <- ggplot(Sum7,                ### The data frame to use.
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr7") +
  
  annotate("text",
           x = 1:3,
           y = Sum3$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum3$.group))

plot8 <- ggplot(Sum8,                ### The data frame to use.
                aes(x = Age,
                    y = emmean)) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05,
                size  = 0.5) +
  geom_point(shape = 15,
             size  = 4) +
  theme_bw() +
  theme(axis.title   = element_text(face  = "bold")) +
  ylab("Proportion of CG methylation Chr8") +
  
  annotate("text",
           x = 1:3,
           y = Sum4$asymp.UCL + 0.008,
           label = gsub(" ", "", Sum4$.group))


