plot2DidLineage <- function( LineageTree, lineageId){

#this snippet requires previously open png with cell plot at plotVals coordinates

lineageCoords 	<- LineageTree[[ lineageId ]]
lineage		<- colnames( lineageCoords)

linesData <- t( lineageCoords)
lines( linesData) 
points( linesData, pch = 16)
}


