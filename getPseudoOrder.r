getPseudoOrder <- function(slingObj){

#this script returns a list of cells

#curveIds 	<- which(slingObj@lineageControl$end.given)
curveIds	<- seq( 1, length( slingObj@lineageControl))
curveNames	<- paste0("curve", curveIds)
curves		<- sapply(curveNames, function(curveName) {keep <- !is.na( psTime[, curveName]);  sort(psTime[ keep, curveName])})

}


