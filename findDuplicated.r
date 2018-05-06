
findDuplicated <- function(CellTable){

source("R/writeDupFile.r")	

	dupInd 		<- which(duplicated(t(CellTable$Genes)))
	dupTable	<- apply(CellTable$Genes[dupInd], 2, function(y) {which(apply(CellTable$Genes,2, function(x) all(x == y)))})

	dupRecord		<- data.frame( 	unlist(CellTable$Cells["FileName",dupTable[1,]]), 
					unlist(CellTable$Cells["batch", dupTable[1,]]),
					unlist(CellTable$Cells["num", dupTable[1,]]),
					unlist(CellTable$Cells["FileName", dupTable[2,]]), 
					unlist(CellTable$Cells["batch", dupTable[2,]]), 
					unlist(CellTable$Cells["num", dupTable[2,]])	)

	names(dupRecord)	<- c("File1", "batch1", "num1", "File2", "batch2", "num2")

	dupFileName		<- paste0( ResPath, .Platform$file.sep, "Duplicates", ".csv" )
	dupFile			<- file(dupFileName, open = "w")
	write.table( dupRecord, file = dupFile, sep = "\t", row.names = FALSE)
	close(dupFile)

	dupTable
}				#main
