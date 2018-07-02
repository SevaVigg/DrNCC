getLineageCoords <- function(seuratObj, slingObj){

cellClust	<- slingObj@clusterLabels	#vector containng cluster ids named with cognate cells
lineageNames	<- names(slingObj@lineages)
lineageCoords	<- sapply( lineageNames, function(x){
 			sapply( slingObj@lineages[[x]], 
				function(y) apply( seuratObj@dr$tsne@cell.embeddings[ names(cellClust)[cellClust == y], ], 2, mean)) })
return(lineageCoords)
}
