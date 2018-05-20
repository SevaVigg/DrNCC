getLineageCoords <- function(seuratObj, slingObj){

lineageIds	<- seq( 1, length( slingObj@lineageControl$end.clus))
lineageNames	<- paste0("Lineage", lineageIds)
lineageCoords	<- sapply( lineageNames, function(x){
 			sapply( slingObj@lineages[[x]], 
				function(y) apply( seuratObj@dr$tsne@cell.embeddings[ WhichCells(seuratObj, as.numeric(y)), ], 2, mean)) })
return(lineageCoords)
}
