getClusterTypes <- function(clFactor){

#This funciton takes a factor whose names are cells and whose values are clusters and assigns colors to clusters with predominant cell types
if(!require(data.table)){
  install.packages("data.table")}

cellClusters <- sapply( unique(levels(clFactor)), function(x) unlist( transpose( strsplit( names( clFactor[clFactor == x]), "\\." ))[[1]]))
nClust	     <- max( as.numeric( levels( clFactor)))


cl_IP	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("IP", x)])==0) 0 else ( 1 + length( grep("IP", x)))/sqrt(1+length(x))))
cl_MC	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("MC", x)])==0) 0 else ( 1 + length( grep("MC", x)))/sqrt(1+length(x))))
cl_tail	<- which.max( sapply(cellClusters, function(x) if(length(x[grep("tails", x)])==0) 0 else ( 1 + length( grep("tails", x)))/sqrt(1+length(x))))

clusterTypes 			<- seq( 1:nClust)
if (length(unique(c(cl_IP, cl_MC, cl_tail))) < 3) {cat( "Some reference clusters conincide \n")
							 attr(clusterType, which = "success") <- TRUE
							 return( clusterTypes)}
names( clusterTypes)		<- rep( "other", nClust)
names( clusterTypes)[ as.numeric(names(cl_IP))]		<- "I"
names( clusterTypes)[ as.numeric(names(cl_MC))]		<- "M"
names( clusterTypes)[ as.numeric(names(cl_tail))]	<- "Tl"

attr(clusterTypes, which = "success") <- TRUE

return(clusterTypes)
}
