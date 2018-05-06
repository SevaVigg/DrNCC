source("R/ReadSourceFiles.r")
source("R/writeDupFile.r")

StartDir <- getwd()

RawPath 	<- file.path( StartDir, "Source", "Raw", "Zerbra_fish_cells")
ResPath		<- file.path( StartDir, "Res")

#Read and process files
CellTable 	<- ReadSourceFiles(RawPath)

write.table(CellTable$Genes, file = "Res/CellTableUnDup_ET.csv", sep = "\t")
write.table(CellTable$Cells, file = "Res/CellTableUnDup_CD.csv", sep = "\t")

#Find duplicates
dupTable  	<- findDuplicated(CellTable)						#Save duplicated entries
writeDupFile(CellTable, dupTable, ResPath)						#my function

			
CellTable$Genes 	<- CellTable$Genes[-dupTable[2,]]
CellTable$Cells 	<- CellTable$Cells[-dupTable[2,]]

#Write deduplicated results

write.csv( CellTable$Genes, file = paste0("Res", .Platform$file.sep, "expressionTableDedup.csv") )
write.csv( CellTable$Cells, file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedup.csv") )
write.csv( CellTable$Probes, file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv") )

genesMissing_I <- which(!complete.cases(CellTable$Genes))
write.table( CellTable$Probes[genesMissing_I,], file = paste0("Res", .Platform$file.sep, "missingGenes.csv") )
lapply(genesMissing_I, function(gene_I) {cell_I <- which(is.na(CellTable$Genes[gene_I,]))
					geneFileName <- paste0("Res", .Platform$file.sep, "Missing_", CellTable$Probes[gene_I, "Gene.Name"], "_cells.csv")
					write.table( CellTable$Cells[,cell_I], file = geneFileName )})




