getClusterTypes <- function(clFactor){

#This funciton takes a factor whose names are cells and whose values are clusters and assigns colors to clusters with predominant cell types
if(!require(data.table)){
  install.packages("data.table")}

cellClusters <- sapply( unique(levels(clFactor)), function(x) unlist( transpose( strsplit( names( clFactor[clFactor == x]), "\\." ))[[1]]))
nClust	     <- length( cellClusters)

cl_IP	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("IP", x)])==0) 0 else ( 1 + length( grep("IP", x)))/sqrt(1+length(x))))
cl_MC	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("MC", x)])==0) 0 else ( 1 + length( grep("MC", x)))/sqrt(1+length(x))))
cl_tail	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("tails", x)])==0) 0 else ( 1 + length( grep("tails", x)))/sqrt(1+length(x))))

clusterTypes 			<- seq( 1:nClust)
names( clusterTypes)		<- rep( "other", nClust)
names( clusterTypes)[ cl_IP] 	<- "I"
names( clusterTypes)[ cl_MC]	<- "M"
names( clusterTypes)[ cl_tail]  <- "Tl"

return(clusterTypes)
}
