

writeDupFile <- function(CellTable, dupTable, ResPath){

dupRecord		<- data.frame( 	CellTable$Cells["FileName",dupTable[1,]], 
					CellTable$Cells["batch", dupTable[1,]],
					CellTable$Cells["num", dupTable[1,]],
					CellTable$Cells["FileName", dupTable[2,]], 
					CellTable$Cells["batch", dupTable[2,]], 
					CellTable$Cells["num", dupTable[2,]]	)

names(dupRecord)	<- c("File1", "batch1", "num1", "File2", "batch2", "num2")

dupFileName		<- paste0( ResPath, .Platform$file.sep, "Duplicates", ".csv" )
dupFile			<- file(dupFileName, open = "w")
write.table( dupRecord, file = dupFile, sep = "\t", row.names = FALSE)
close(dupFile)

} #main
