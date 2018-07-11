getLineageCoords <- function(seuratObj, slingObj){

cellClust	<-  sapply(rownames(slingObj@clusterLabels), function(x) colnames(slingObj@clusterLabels)[which(as.logical(slingObj@clusterLabels[x, ]))])
		#vector containng cluster ids named with cognate cells

lineageNames	<- names(slingObj@lineages)
lineageCoords	<- sapply( lineageNames, function(x){
 			sapply( slingObj@lineages[[x]], 
				function(y) apply( seuratObj@dr$tsne@cell.embeddings[ names(cellClust)[cellClust == y], ], 2, mean)) })
return(lineageCoords)
}
