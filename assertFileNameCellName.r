assertFileNameCellName <- function( cell_list, FileProperties){

assertVar <- function(Var) {	ans_attr 	<- list(1) 
				names(ans_attr) <- "type"
				attributes(assertVar) <- ans_attr

				if( is.na( FileProperties[[Var]]) & !is.na( cell_list[[Var]])){ 
						if( Var != "num"){cat("In FileName ", Var, "is missing, updated from CellName ", "\n")}
						assertVar <- NA; attr(assertVar, "type") <-"take_Cell"
						return(assertVar) }else{
				if( is.na( cell_list[[Var]]) & !is.na( FileProperties[[Var]])){
						cat("In CellName ", Var, "is missing, updated from FileName", "\n" ) 
						assertVar <- NA; attr(assertVar, "type") <- "take_File"
						return(assertVar) }else{
				if( is.na( cell_list[[Var]]) & is.na( FileProperties[[Var]])){
						if( Var != "dateEx"){
						   cat( Var, " is missing both in CellName and FileName", "\n" )
						   assertVar <- FALSE; attr( assertVar, "type") <- paste("Missing CellName and FileName ", Var," values") 
						   return(assertVar)
						}else{
						   cat( Var, "is missing both in CellName and FileName", "\n")
						   assertVar <- NA; attr(assertVar, "type") <-"take_zero"
						   return(assertVar) 	
						} 
				}else{
				if(FileProperties[[Var]] != cell_list[[Var]]) {
						cat( "Warning. CellName ", cell_list[[Var]],  " has ", Var," ", cell_list[[Var]], 
						     " which does not correspond to FileName's ", Var, " ", FileProperties[[Var]], " ")
						if( Var == "CellType"){ assertVar <- NA; attr( assertVar, "type") <- "take_File"; cat( " File type is used\n")}else{
							assertVar <- NA; attr( assertVar, "type") <- "take_Cell"; cat(" Cell type is used\n")}
						return(assertVar) 
									}else{ return(TRUE)}}}}
							} 

assertRes	<- lapply( names(FileProperties), function(x) assertVar(x))
names(assertRes)<- names(FileProperties)

assertRes
}

