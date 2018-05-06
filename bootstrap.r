wDir	<-	getwd()
btsTbl	<-	read.table( "Source/Bootstrap/again/bstCl_21.txt", row.names = 1 , sep = " " )

Cells	<- 	read.table( "Source/Bootstrap/again/CellClustersAll3.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

prmLst	<-	scan( "Source/Bootstrap/again/table_permutation.txt", what = "", sep = "\n" )
prmTbl	<-	strsplit( prmLst, "[[:space:]]+")
prmTbl	<-	lapply(prmTbl, as.integer)
 
