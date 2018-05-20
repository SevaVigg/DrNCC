getLineageTree <- function(slingObj){

lineageIds	<- seq( 1, length( slingObj@lineageControl$end.clus))
lineageNames	<- paste0("Lineage", lineageIds)
lineages	<- sapply(lineageNames, function(x) slingObj@lineages[[x]])
return(lineages)
}
