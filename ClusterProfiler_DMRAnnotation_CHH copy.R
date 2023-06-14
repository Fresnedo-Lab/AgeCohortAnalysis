library(clusterProfiler)
library(rols)
library(ggplot2)
BiocManager::install("DOSE")
library(DOSE)

setwd("~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/")
GoMap <- read.table("GOmap.txt", header = TRUE)
s <- strsplit(GoMap$GOterm, split = ",")
GoMap.new <- data.frame(GeneID = rep(GoMap$GeneID, sapply(s, length)), GOterm = unlist(s))
GoMap.new <- GoMap.new[,c("GOterm", "GeneID")]


GeneList0812 <- read.table("Heatmap/VennDiagrams/CHH/Hyper/CHH_0812_geneList_Hyp.txt", header = TRUE)
GeneList0817 <- read.table("Heatmap/VennDiagrams/CHH/Hyper/CHH_0817_geneList_Hyp.txt", header = TRUE)
GeneList1217 <- read.table("Heatmap/VennDiagrams/CHH/Hyper/CHH_1217_geneList_Hyp.txt", header = TRUE)

DMR_GO_0812 <- enricher(GeneList0812$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)
DMR_GO_0817 <- enricher(GeneList0817$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)
DMR_GO_1217 <- enricher(GeneList1217$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)


head(DMR_GO_0812)
head(DMR_GO_0817)
head(DMR_GO_1217)


head(DMR_GO_0812@result)
head(DMR_GO_0817@result)
head(DMR_GO_1217@result)

GeneList0812_Hypo <- read.table("Heatmap/VennDiagrams/CHH/Hypo/CHH_0812_geneList_Hypo.txt", header = TRUE)
GeneList0817_Hypo <- read.table("Heatmap/VennDiagrams/CHH/Hypo/CHH_0817_geneList_Hypo.txt", header = TRUE)
GeneList1217_Hypo <- read.table("Heatmap/VennDiagrams/CHH/Hypo/CHH_1217_geneList_Hypo.txt", header = TRUE)

DMR_GO_0812_Hypo <- enricher(GeneList0812_Hypo$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)
DMR_GO_0817_Hypo <- enricher(GeneList0817_Hypo$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)
DMR_GO_1217_Hypo <- enricher(GeneList1217_Hypo$GeneID, pvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = GoMap.new)


head(DMR_GO_0812_Hypo)
head(DMR_GO_0817_Hypo)
head(DMR_GO_1217_Hypo)


head(DMR_GO_0812_Hypo@result)
head(DMR_GO_0817_Hypo@result)
head(DMR_GO_1217_Hypo@result)

Merged_GO_0812 <- merge(GeneList0812, GoMap.new, by = "GeneID")
Merged_GO_0817 <- merge(GeneList0817, GoMap.new, by = "GeneID")
Merged_GO_1217 <- merge(GeneList1217, GoMap.new, by = "GeneID")

Merged_GO_0812_Hypo <- merge(GeneList0812_Hypo, GoMap.new, by = "GeneID")
Merged_GO_0817_Hypo <- merge(GeneList0817_Hypo, GoMap.new, by = "GeneID")
Merged_GO_1217_Hypo <- merge(GeneList1217_Hypo, GoMap.new, by = "GeneID")

go <- Ontology("go")

df <- data.frame()
for (i in 1:nrow(Merged_GO_0812)){
  GO <- Merged_GO_0812$GOterm[i]
  geneID <- Merged_GO_0812$GeneID[i]
  term <- term(go, GO)
  label <- termLabel(term)
  annot <- as.character(term@annotation$has_obo_namespace[1])
  entry <- as.data.frame(cbind(geneID, GO, annot, label))
  df <- rbind(df, entry)
}

#CHH_Hyper
df_biopro_0812 <- df[ which(df$annot=='biological_process'), ]
df_molfun_0812 <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_0812 <- df[ which(df$annot=='cellular_component'), ]

bioproc_0812 <- aggregate(data.frame(count = df_biopro_0812$label), list(value = df_biopro_0812$label), length)
CHH_bioproc_0812 <- bioproc_0812[order(-bioproc_0812$count),]
write.csv(CHH_bioproc_0812, "CHH_BioProc_0812_hyper.csv")
molfun_0812 <- aggregate(data.frame(count = df_molfun_0812$label), list(value = df_molfun_0812$label), length)
CHH_molfun_0812 <- molfun_0812[order(-molfun_0812$count),]
write.csv(CHH_molfun_0812, "CHH_MolFun_0812_hyper.csv")
cellcom_0812 <- aggregate(data.frame(count = df_cellcom_0812$label), list(value = df_cellcom_0812$label), length)
CHH_cellcom_0812 <- cellcom_0812[order(-cellcom_0812$count),]
write.csv(CHH_cellcom_0812, "CHH_CellCom_0812_hyper.csv")

df_biopro_0817 <- df[ which(df$annot=='biological_process'), ]
df_molfun_0817 <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_0817 <- df[ which(df$annot=='cellular_component'), ]

bioproc_0817 <- aggregate(data.frame(count = df_biopro_0817$label), list(value = df_biopro_0817$label), length)
CHH_bioproc_0817 <- bioproc_0817[order(-bioproc_0817$count),]
write.csv(CHH_bioproc_0817, "CHH_BioProc_0817_hyper.csv")
molfun_0817 <- aggregate(data.frame(count = df_molfun_0817$label), list(value = df_molfun_0817$label), length)
CHH_molfun_0817 <- molfun_0817[order(-molfun_0817$count),]
write.csv(CHH_molfun_0817, "CHH_MolFun_0817_hyper.csv")
cellcom_0817 <- aggregate(data.frame(count = df_cellcom_0817$label), list(value = df_cellcom_0817$label), length)
CHH_cellcom_0817 <- cellcom_0817[order(-cellcom_0817$count),]
write.csv(CHH_cellcom_0817, "CHH_CellCom_0817_hyper.csv")

df_biopro_1217 <- df[ which(df$annot=='biological_process'), ]
df_molfun_1217 <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_1217 <- df[ which(df$annot=='cellular_component'), ]

bioproc_1217 <- aggregate(data.frame(count = df_biopro_1217$label), list(value = df_biopro_1217$label), length)
CHH_bioproc_1217 <- bioproc_1217[order(-bioproc_1217$count),]
write.csv(CHH_bioproc_1217, "CHh_BioProc_1217_hyper.csv")
molfun_1217 <- aggregate(data.frame(count = df_molfun_1217$label), list(value = df_molfun_1217$label), length)
CHH_molfun_1217 <- molfun_1217[order(-molfun_1217$count),]
write.csv(CHH_molfun_1217, "CHH_MolFun_1217_hyper.csv")
cellcom_1217 <- aggregate(data.frame(count = df_cellcom_1217$label), list(value = df_cellcom_1217$label), length)
CHH_cellcom_1217 <- cellcom_1217[order(-cellcom_1217$count),]
write.csv(CHH_cellcom_1217, "CHH_CellCom_1217_hyper.csv")

#CHG_Hypo
df_biopro_0812_hypo <- df[ which(df$annot=='biological_process'), ]
df_molfun_0812_hypo <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_0812_hypo <- df[ which(df$annot=='cellular_component'), ]

bioproc_0812_hypo <- aggregate(data.frame(count = df_biopro_0812_hypo$label), list(value = df_biopro_0812_hypo$label), length)
CHH_bioproc_0812_hypo <- bioproc_0812_hypo[order(-bioproc_0812_hypo$count),]
write.csv(CHH_bioproc_0812_hypo, "CHH_BioProc_0812_hypo.csv")
molfun_0812_hypo <- aggregate(data.frame(count = df_molfun_0812_hypo$label), list(value = df_molfun_0812_hypo$label), length)
CHH_molfun_0812_hypo <- molfun_0812_hypo[order(-molfun_0812_hypo$count),]
write.csv(CHH_molfun_0812_hypo, "CHH_MolFun_0812_hypo.csv")
cellcom_0812_hypo <- aggregate(data.frame(count = df_cellcom_0812_hypo$label), list(value = df_cellcom_0812_hypo$label), length)
CHH_cellcom_0812_hypo <- cellcom_0812_hypo[order(-cellcom_0812_hypo$count),]
write.csv(CHH_cellcom_0812_hypo, "CHH_CellCom_0812_hypo.csv")

df_biopro_0817_hypo <- df[ which(df$annot=='biological_process'), ]
df_molfun_0817_hypo <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_0817_hypo <- df[ which(df$annot=='cellular_component'), ]

bioproc_0817_hypo <- aggregate(data.frame(count = df_biopro_0817_hypo$label), list(value = df_biopro_0817_hypo$label), length)
CHH_bioproc_0817_hypo <- bioproc_0817_hypo[order(-bioproc_0817_hypo$count),]
write.csv(CHH_bioproc_0817_hypo, "CHH_BioProc_0817_hypo.csv")
molfun_0817_hypo <- aggregate(data.frame(count = df_molfun_0817_hypo$label), list(value = df_molfun_0817_hypo$label), length)
CHH_molfun_0817_hypo <- molfun_0817_hypo[order(-molfun_0817_hypo$count),]
write.csv(CHH_molfun_0817_hypo, "CHH_MolFun_0817_hypo.csv")
cellcom_0817_hypo <- aggregate(data.frame(count = df_cellcom_0817_hypo$label), list(value = df_cellcom_0817_hypo$label), length)
CHH_cellcom_0817_hypo <- cellcom_0817_hypo[order(-cellcom_0817_hypo$count),]
write.csv(CHH_cellcom_0817_hypo, "CHH_CellCom_0817_hypo.csv")

df_biopro_1217_hypo <- df[ which(df$annot=='biological_process'), ]
df_molfun_1217_hypo <- df[ which(df$annot=='molecular_function'), ]
df_cellcom_1217_hypo <- df[ which(df$annot=='cellular_component'), ]

bioproc_1217_hypo <- aggregate(data.frame(count = df_biopro_1217_hypo$label), list(value = df_biopro_1217_hypo$label), length)
CHH_bioproc_1217_hypo <- bioproc_1217_hypo[order(-bioproc_1217_hypo$count),]
write.csv(CHH_bioproc_1217_hypo, "CHH_BioProc_1217_hypo.csv")
molfun_1217_hypo <- aggregate(data.frame(count = df_molfun_1217_hypo$label), list(value = df_molfun_1217_hypo$label), length)
CHH_molfun_1217_hypo <- molfun_1217_hypo[order(-molfun_1217_hypo$count),]
write.csv(CHH_molfun_1217_hypo, "CHH_MolFun_1217_hypo.csv")
cellcom_1217_hypo <- aggregate(data.frame(count = df_cellcom_1217_hypo$label), list(value = df_cellcom_1217_hypo$label), length)
CHH_cellcom_1217_hypo <- cellcom_1217_hypo[order(-cellcom_1217_hypo$count),]
write.csv(CHH_cellcom_1217_hypo, "CHH_CellCom_1217_hypo.csv")
