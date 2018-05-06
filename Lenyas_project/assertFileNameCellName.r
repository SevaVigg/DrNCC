assertFileNameCellName <- function( cell_lists, FileProperties){
#It's must check hpf field in csv
ans <- TRUE
	lapply( cell_lists, function(x) { 
		if(is.na(FileProperties$hour) || (FileProperties$hour != x$hour))
		{ 
			cat( "Error. Cell ", x$num,  " has hpf ", x$hour, 
			" which does not correspond to file hpf ", FileProperties$hour, "\n") 
		} } )

ans
}

