input.dir<- "~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/"
output.dir<- "~/Documents/Almond/AgeCohortProject/DSSAnalysis_DMRCalling/"

#read data
setwd(input.dir)
dat<- read.table("Almond_Feb_2021_with_functional_info.gff", comment.char="#", sep="\t")
#Keep only RNA rows
mRNA<- subset(dat, V3=="mRNA")
#Keep only rows that include GO terms
mRNA.go<- dplyr::filter(mRNA, grepl("Ontology_term",V9))
ids<- as.vector(mRNA.go$V9)
#Isolate GO terms using regular expressions
GO<- sub(".+Ontology_term=([GO:,0123456789]+).*", "\\1", ids)
#Isolate mRNA names using regular expressions
mRNA_name<- sub(".+Name=", "\\1", ids)
mRNA_name2 <- sub(";.+", "", gene)
#Bonus: If desired, you can also isolate the InterPro IDs
ipro<- sub(".+InterPro:(IPR[0123456789]+).+", "\\1", ids)

#Make map by combining GO term a gene names
GOmap<- data.frame(cbind(GO, mRNA_name2), stringsAsFactors=F)
GOmap <- GOmap[,c("mRNA_name2", "GO")]
#Make one GO term-by-gene per row
GOlist<- tidyr::separate_rows(GOmap, "GO", sep=",")

#Write file to use as GO2GENE in clusterProfiler::enricher
setwd(output.dir)
write.table(GOlist, "GOlist.txt", row.names=F, quote=F, sep="\t")
write.table(GOmap, "GOmap.txt", row.name = F, quote = F, sep="\t")
