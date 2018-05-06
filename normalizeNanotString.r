require(NanoStringNorm)
require(preprocessCore)



Genes 	<- read.csv( file = paste0("Res", .Platform$file.sep, "expressionTableDedupQC.csv") )
Cells	<- read.csv( file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedupQC.csv") )
Probes	<- read.csv( file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv") )

geneNames <- Genes[,1]
cellNames <- Cells[,1]

Cells	<- as.data.frame(lapply(Cells, as.character), stringsAsFactors = FALSE)
Genes	<- as.data.frame(lapply(Genes, as.integer), stringsAsFactors = FALSE)
Probes	<- as.data.frame(lapply(Probes, as.character), stringsAsFactors = FALSE)

rownames(Cells)	<- cellNames
rownames(Genes)	<- geneNames
Cells[,1]	<- NULL
Genes[,1]	<- NULL

NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes, stringsAsFactors = FALSE)
colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")

NanoTable[c("Kanamycin Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
NanoTableNormed <- NanoStringNorm(x = NanoTable, CodeCount = "sum", Background = "mean", SampleContent = "housekeeping.geo.mean")



