plot2DallLineages <- function( LineageTree, plotVals, clustTypes, plotDir, LineageId){

source("R/setClusterColor.r") 
cellColors <- setClusterColor(ipmc, clustTypes, LineageTree[[LineageId]])

png( file.path( plotDir, "allLineages.png"))
	plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
	sapply(LineageTree, function(x) {lines(t(x)); points(t(x), pch = 16)})
dev.off()
}


