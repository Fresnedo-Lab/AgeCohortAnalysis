BiocManager::install("topGO")
BiocManager::install("ALL")
BiocManager::install("Rgraphviz")
library(topGO)

setwd("~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/")

geneID2GO_Alm <- readMappings("GOmap2.txt")
head(geneID2GO_Alm)
str(geneID2GO_Alm)
geneNames_Alm <- names(geneID2GO_Alm)
head(geneNames_Alm)
CG_0812_HypGenes <- read.delim("Heatmap/VennDiagrams/CG_0812_geneList_Hyp.txt", header = FALSE)
CG_0812_HypGenes <- as.vector(CG_0812_HypGenes$V1)
geneList_Alm <- factor(as.integer(geneNames_Alm %in% CG_0812_HypGenes))

CG_0817_HypGenes <- read.delim("Heatmap/VennDiagrams/CG_0817_geneList_Hyp.txt", header = FALSE)
CG_0817_HypGenes <- as.vector(CG_0817_HypGenes$V1)
geneList_Alm0817 <- factor(as.integer(geneNames_Alm %in% CG_0817_HypGenes))

CG_1217_HypGenes <- read.delim("Heatmap/VennDiagrams/CG_1217_geneList_Hyp.txt", header = FALSE)
CG_1217_HypGenes <- as.vector(CG_1217_HypGenes$V1)
geneList_Alm1217 <- factor(as.integer(geneNames_Alm %in% CG_1217_HypGenes))

names(geneList_Alm) <- geneNames_Alm
str(geneList_Alm)

names(geneList_Alm0817) <- geneNames_Alm
str(geneList_Alm0817)

names(geneList_Alm1217) <- geneNames_Alm
str(geneList_Alm1217)

GOdata_CG_0812_Hyp <- new("topGOdata", ontology = "BP", allGenes = geneList_Alm, annot = annFUN.gene2GO, gene2GO = geneID2GO_Alm)
GOdata_CG_0817_Hyp <- new("topGOdata", ontology = "BP", allGenes = geneList_Alm0817, annot = annFUN.gene2GO, gene2GO = geneID2GO_Alm)
GOdata_CG_1217_Hyp <- new("topGOdata", ontology = "BP", allGenes = geneList_Alm1217, annot = annFUN.gene2GO, gene2GO = geneID2GO_Alm)

resultsFisher <- runTest(GOdata_CG_0812_Hyp, algorithm = "classic", statistic = "fisher")
resultsFisher
resultKS <- runTest(GOdata_CG_0812_Hyp, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata_CG_0812_Hyp, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata_CG_0812_Hyp, classicFisher = resultsFisher, classicKS = resultKS, elimKS = resultKS.elim,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes
write.table(allRes, "CG_Hyp_0812.txt")
showSigOfNodes(GOdata_CG_0812_Hyp, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

resultsFisher0817 <- runTest(GOdata_CG_0817_Hyp, algorithm = "classic", statistic = "fisher")
resultsFisher0817
resultKS0817 <- runTest(GOdata_CG_0817_Hyp, algorithm = "classic", statistic = "ks")
resultKS.elim0817 <- runTest(GOdata_CG_0817_Hyp, algorithm = "elim", statistic = "ks")

allRes0817 <- GenTable(GOdata_CG_0817_Hyp, classicFisher = resultsFisher0817, classicKS = resultKS0817, elimKS = resultKS.elim0817,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes0817
write.table(allRes0817, "CG_Hyp_0817.txt")
showSigOfNodes(GOdata_CG_0817_Hyp, score(resultKS.elim0817), firstSigNodes = 5, useInfo = 'all')

resultsFisher1217 <- runTest(GOdata_CG_1217_Hyp, algorithm = "classic", statistic = "fisher")
resultsFisher1217
resultKS1217 <- runTest(GOdata_CG_1217_Hyp, algorithm = "classic", statistic = "ks")
resultKS.elim1217 <- runTest(GOdata_CG_1217_Hyp, algorithm = "elim", statistic = "ks")

allRes1217 <- GenTable(GOdata_CG_1217_Hyp, classicFisher = resultsFisher1217, classicKS = resultKS1217, elimKS = resultKS.elim1217,orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)
allRes1217
write.table(allRes1217, "CG_Hyp_1217.txt")
showSigOfNodes(GOdata_CG_1217_Hyp, score(resultKS.elim1217), firstSigNodes = 5, useInfo = 'all')


