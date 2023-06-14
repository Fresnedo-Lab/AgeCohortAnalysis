input.dir<- "~/Documents/CAPS/RNAseq/Tomato/genome/Heinz1706/annotation/ITAG4.1_release"
output.dir<- "~/Documents/CAPS/RNAseq/Tomato/"

#read data
setwd(input.dir)
dat<- read.table("ITAG4.1_gene_models.gff", comment.char="#", sep="\t")
#Keep only RNA rows
RNAs<- subset(dat, V3=="mRNA")
#Keep only rows that include GO terms
RNAs.go<- dplyr::filter(RNAs, grepl("Ontology_term",V9))
ids<- as.vector(RNAs.go$V9)
#Isolate GO terms using regular expressions
GO<- sub(".+Ontology_term=([GO:,0123456789]+).*", "\\1", ids)
#Isolate mRNA names using regular expressions
mRNA<- sub(".+Name=(Solyc[g.0123456789]+).+", "\\1", ids)
#Bonus: If desired, you can also isolate the InterPro IDs
ipro<- sub(".+InterPro:(IPR[0123456789]+).+", "\\1", ids)
#Gene name is mRNA name minus last two characters
gene<- substr(mRNA, 1, nchar(mRNA)-2)

#Make map by combining GO term a gene names
GOmap<- data.frame(cbind(GO, gene), stringsAsFactors=F)
#Make one GO term-by-gene per row
GOlist<- tidyr::separate_rows(GOmap, "GO", sep=",")

#Write file to use as GO2GENE in clusterProfiler::enricher
setwd(output.dir)
write.table(GOlist, "GOlist.txt", row.names=F, quote=F, sep="\t")
