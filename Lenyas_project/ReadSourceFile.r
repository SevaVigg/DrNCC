source("ParseNanoStringName.r")
source("assertFileNameCellName.r")
source("classNames.r")
library(plyr)


getTables <- function (filename, FileProperties) {
	datafile 	<- file( filename, open = "r" ) 
	data 		<- read.csv( datafile, sep=",", header=TRUE, stringsAsFactors=FALSE, check.names = FALSE )
	close(datafile)

	cell_inds	<- grep( "study", colnames( data ) )	
	if (length(cell_inds) == 0) {
		cell_inds	<- grep( "hpf", colnames( data ) )
	}
	cell_strs 	<- colnames( data)[cell_inds]
	cell_lists	<- lapply( cell_strs, ParseNanoStringName)
	#We try to fix columns from filename (if it's possible)
	cell_lists  <- lapply(cell_lists, function(x) {
			if (is.na(x$hour)) {
				x$hour <- FileProperties$hour
			}
			return(x)
		}
	)
	if( !assertFileNameCellName( cell_lists, FileProperties)){ return(ans <- NA)}
	data <- data[-1,]
	Descs <- as.data.frame( do.call( cbind, cell_lists))
	resTable <- data[, cell_inds]
	#We have 61 genes and controls. If it's too much, we must warn user.  
	if (length(resTable[,1]) != 61) {
		cat("File ", filename, "haven't all nesessary genes or have extra genes!\n")
		return(NA)
	}
	colnames(resTable) <- paste(unlist( Descs["date",]),".", unlist( Descs["hour",]),".",unlist( Descs["num",]), sep = "")
	colnames(Descs)	   <- colnames(resTable)
	#Stop if we don't have GeneName column
	if (!("Gene Name" %in% colnames(data))) {
		return(NA)
	}
	rownames(resTable) <- data[, "Gene Name"]
	if (!("Target Sequence" %in% colnames(data))) {
		data[,"Target Sequence"] <- sprintf("",seq(1:length(rownames(data))))
	}
	if (!("Annotation" %in% colnames(data))) {
		data[,"Annotation"] <- sprintf("",seq(1:length(rownames(data))))
	}
	if (!("Class Name" %in% colnames(data))) {
		classNameCurrSet <- laply(data[,"Gene Name"], function (x) {
				return(className[[x]])
			}
		)
		data[,"Class Name"] <- classNameCurrSet
	}
	GenesProbes	<- data[,c("Gene Name", "Annotation", "Accession #", "Class Name", "Target Sequence")]

	ans		<- list()	 	
	ans$genes	<- resTable	
	ans$Descs	<- Descs
	ans$GenesProbes	<- GenesProbes	
	ans
}

parseFile <- function (FileName) {
	FileProperties <- ParseNanoStringName(tail(strsplit(FileName, "/")[[1]], n=1))
	resTable       <- getTables(FileName, FileProperties )	
	return(resTable)
}
#!!IMPORTANT!!
#Here is csv-files directory. You must change it, when you move this project to another directory
StartDir <- "/home/leo/fish/data"
SourcePath <- dir( StartDir )

filenames <- list.files(path = paste(StartDir, SourcePath, sep = "/") , pattern=".csv$", full.names = TRUE )
#Open all files in directory
tables <- lapply(filenames, parseFile)
tables<-tables[!is.na(tables)]
library(purrr)
geneTable<-lapply(tables, function(x) {
	return(as.data.frame(x$genes))
})
#Collect the summary table from list of tables
unitedTable$genes <- reduce(geneTable, cbind)

