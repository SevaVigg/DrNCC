assignColors <- function(seurObj){

#This funciton takes a seurat object and assigns colors to clusters with predominant cell types

source("R/getClusterTypes.r")

clusterTypes <- getClusterTypes(seurObj)

cellColors <- rep("grey", length( clusterTypes))

cellColors[ clusterTypes["M"]]	<-   	"black"
cellColors[ clusterTypes["Tl"]]	<- 	"red"
cellColors[ clusterTypes["I"]]	<-   	"cyan"

return(cellColors)
}
