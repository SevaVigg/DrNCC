getPseudoOrder <- function(slingObj){

#this script returns a list of cells

#curveIds 	<- which(slingObj@lineageControl$end.given)
curveIds	<- seq_along( slingObj@lineageControl$end.clus)
curveNames	<- paste0("curve", curveIds)
curves		<- sapply(curveNames, function(curveName) {keep <- !is.na( psTime[, curveName]);  sort(psTime[ keep, curveName])})

}


