ReadSourceFiles <- function(SourcePath){

# Genes containing gene descriptions. 

source("R/ReadNanoStringFile.r")
source("R/findDuplicated.r")
source("R/writeDupFile.r")

RawPath		<- file.path( SourcePath, "All")
#FilteredPath	<- file.path( SourcePath, "Filtered")
BadPath		<- file.path( SourcePath, "BadFiles")
ResPath		<- file.path( SourcePath, "Res")

dir.create( BadPath, showWarnings = FALSE)
dir.create( ResPath, showWarnings = FALSE)

filenames <- list.files(path = RawPath, pattern=".csv$", full.names = TRUE )

for( FileName in filenames){								#find the first correct file
	CellTable       <- ReadNanoStringFile( FileName)
	filenames	<- filenames[-1]	
	if(is.na( CellTable)){
		file.rename( from = FileName,  to = file.path( BadPath, basename(FileName))); 
		next }else{ break}
}


for( FileName in filenames){
	NewTable	<- ReadNanoStringFile( FileName)
	filenames	<- filenames[-1]
	if(is.na(NewTable)){file.rename(from = FileName,  to = file.path( BadPath, basename(FileName))); next }

	
	if( !all(unlist(lapply( rownames(NewTable$Cells), function(Var) {
			if( all(!is.element(NewTable$Cells[Var,], CellTable$Cells[Var,])) 
			    & identical(CellTable$Genes, NewTable$Genes)){
			return( TRUE)}else{return(FALSE)}
								})))){
			CellTable$Genes <- cbind(CellTable$Genes, NewTable$Genes)
			CellTable$Cells <- cbind(CellTable$Cells, NewTable$Cells)
	    }else{
		cat(FileName, " already present", "\n")
	}
}

DescNames		<- rownames(CellTable$Cells)

Genes 			<- as.data.frame(lapply(CellTable$Genes, unlist), stringsAsFactors = FALSE)
rownames(Genes)		<- CellTable$Probes[, "Gene Name"] 
CellTable$Cells 	<- as.data.frame(lapply(CellTable$Cells, unlist), stringsAsFactors = FALSE)
CellTable$Cells		<- as.data.frame(lapply(CellTable$Cells, as.character), stringsAsFactors = FALSE)
CellTable$Probes	<- as.data.frame(lapply(CellTable$Probes, as.character), stringsAsFactor = FALSE)
rownames(CellTable$Cells)	<- DescNames

ans <- CellTable

} #main



