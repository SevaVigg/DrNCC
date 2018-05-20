plot2DidLineage <- function( seuratObj, LineageTee, lineageId, MC_linId, IP_linId, plotVals, clustTypes, plotDir){

lineageCoords 	<- LineageTree[[ lineageId ]]
lineage		<- colnames( lineageCoords)
source("R/setLineageColors.r")

cellColors <- rep("gray", length( seuratObj@ident))
cellColors[ which(seuratObj@ident == clustTypes[ "I"]) ]  <- "cyan"
cellColors[ which(seuratObj@ident == clustTypes[ "M"]) ]  <- "black"
lineageColors <-  setLineageColors( lineage, lineageId, MC_linId, IP_linId)

for (lind in seq(1, length(lineage))) cellColors[ which( seuratObj@ident == lineage[lind])] <- lineageColors[lind]
png( file.path( plotDir, paste0("lineage", lineageId, "_plot.png")))
	plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
	linesData <- t( lineageCoords)
	lines( linesData) 
	points( linesData, pch = 16)
dev.off()
}


