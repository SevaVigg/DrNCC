getClusterTypes <- function(seurObj){

#This funciton takes a :seurat object and assigns colors to clusters with predominant cell types

cellClusters <- sapply( unique(levels(seurObj@ident)), function(x) table(unlist(Cells["hpf_CellType" , WhichCells(seurObj, x)])))
nClust	     <- length( cellClusters)

cl_IP	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("IP", names(x))])==0) 0 else x[grep("IP", names(x))]))
cl_MC	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("MC", names(x))])==0) 0 else x[grep("MC", names(x))]))
cl_tail	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("tails", names(x))])==0) 0 else x[grep("tails", names(x))]))

clusterTypes 			<- seq( 1:nClust)
names( clusterTypes)		<- rep( "other", nClust)
names( clusterTypes)[ cl_IP] 	<- "I"
names( clusterTypes)[ cl_MC]	<- "M"
names( clusterTypes)[ cl_tail]  <- "Tl"

return(clusterTypes)
}
