ReadNanoStringFile <- function( FileName){

source("R/assertFileNameCellName.r")
source("R/ParseNanoStringName.r")
source("R/classNames.r")					#the array with names equeal to gene names
#library(plyr)

	FileProperties <- ParseNanoStringName( FileName)

	cat( "Processing file ", FileName, "\n")

	datafile 	<- file( FileName, open = "r" ) 
	data 		<- read.csv( datafile, sep=",", header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
	close(datafile)

	ans		<- list()	 	
	ans_attr 	<- list(1) 
	names(ans_attr) <- "diagn"
	attributes(ans) <- ans_attr	

	cell_inds	<- grep( "RCC", colnames( data ) )

	if( length(cell_inds)==0){ cat( "No columns marked as RCC", "\n"); ans <- NA; attr(ans, "diagn") <- paste( FileName, "has no RCC columns");
				  return(ans)}
		
	cell_strs 	<- colnames( data)[cell_inds]
	cell_lists	<- lapply( cell_strs, ParseNanoStringName)	#ParseNanoStringName creates the dictionary


	#Check column names
	
	#GeneName column stop
	if (!("Gene Name" %in% colnames(data))) { cat( "No gene name column", "\n"); ans <- NA; attr(ans, "diagn") <- paste( FileName, "has no gene name columns")
							return(ans)}

	#Some file contains empty and senseless lines in the beginning. Removing all of them till "Kanamycin Pos"

	Kanamycin_strInd <- which( data[,"Gene Name"]=="Kanamycin Pos")
	ForwardStrInds	 <- if (Kanamycin_strInd == 1) numeric() else seq(1, Kanamycin_strInd -1)	#because seq(1,0) = c(1,0)

	if( length(ForwardStrInds) >0  ){ data <- data[-ForwardStrInds,]
		 cat( "Removing ", length(ForwardStrInds), " empty string in the beginning", FileName, "\n")
		}else{ cat( "No initial strings to remove", "\n")}

	#Now check if all Genes are here

	if( length( setdiff( names(className), data$"Gene Name") )>0 ){
		cat("Genes ", setdiff(names(className), data$"Gene Name"), " are missing in ", FileName, "\n")
		ans <- NA; attr(ans, "diagn") <- paste( FileName, "Incomplete gene set"); return(ans)}
	if( length( setdiff( data$"Gene Name", names(className) ))>0 ){
		cat("Warning: Genes ", setdiff(data$"Gene Name", names(className)), " are redundant in ", FileName, "\n")}

	#remove extra strings, keep only the Genes that are in the className

	data <- data[data$"Gene Name" %in% names(className),]
	
	#If we dont have less important columns try to update

	if (!("Target Sequence" %in% colnames(data))) {
		data[,"Target Sequence"] <- sprintf("",seq(1:nrow(data)))
	}
	if (!("Annotation" %in% colnames(data))) {
		data[,"Annotation"] <- sprintf("",seq(1:nrow(data)))
	}
	if (!("Class Name" %in% colnames(data))) {
		classNameCurrSet <- lapply(data[,"Gene Name"], function (x) {return(className[[x]]) })
		data$"Class Name" <- classNameCurrSet
		cat("Class Name is missing in ", FileName, " updated from className", "\n")
	}

	Cells	<- lapply( cell_lists, function(testCell){
					cat( testCell$"num", "\t")
					asrtCell <- assertFileNameCellName(testCell, FileProperties)
					testCell_ans <- lapply( names(asrtCell), function(Var){
							if( is.na( asrtCell[[Var]] ) ){
								if( attr(asrtCell[[Var]], "type") == "take_zero" ){cat(Var, "updated with zero \n"); return(0)}
								if( attr(asrtCell[[Var]], "type") == "take_Cell" ){
                                                             		return( testCell[[Var]])
                                                           	}else{ return( FileProperties[[Var]])
								}
							}else{ 	if( !asrtCell[[Var]] ){cat("Wrong file ", FileName, "\n"); return( asrtCell[[Var]])
						       		}else{ return( testCell[[Var]])
								}	
							}
											}
						)
					if( length( which( testCell_ans %in% FALSE)) != 0 ){
						return( testCell_ans[[ which( testCell_ans %in% FALSE)[1] ]] )
					}else{return(testCell_ans)}						
								}
				)

	if( length( which( Cells %in% FALSE)) != 0){
		ans <- NA
		attr(ans, "diagn") <- attr(Cells[[ which( Cells %in% FALSE)[1] ]], "type")
		return(ans)
		}
	
	CellDesc	   	<- as.data.frame( do.call(cbind, Cells))
	rownames(CellDesc) 	<- names(FileProperties)
	CellDesc		<- rbind(CellDesc, basename(FileName))	
	rownames(CellDesc)[nrow(CellDesc)] <- "FileName"


	resTable 			<- data[, cell_inds]
	resTable			<- as.data.frame(lapply(resTable, as.numeric))
	rownames(resTable) 		<- names(className)
	colnames(resTable) 		<- unlist( CellDesc["hpf",])
	colnames(CellDesc)	  	<- colnames(resTable)



	Probes	<- data[,c("Gene Name", "Annotation", "Accession #", "Class Name", "Target Sequence")]


	ans$Genes	<- resTable	
	ans$Cells	<- CellDesc
	ans$Probes	<- Probes
	
	ans
}				#main

