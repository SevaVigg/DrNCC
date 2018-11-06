plot2DIPMCLineages <- function( LineageTree, plotVals, clustTypes, plotDir){

source("R/setClusterColors.r")
cellColors <- setClusterColors(ipmc, clustTypes)

endCluster <- sapply(LineageTree, function(x) colnames(x)[length(colnames(x))] )
MC_linId   <- names(endCluster[endCluster==clustTypes["M"]]) 
IP_linId   <- names(endCluster[endCluster==clustTypes["I"]]) 

png( file.path( plotDir, "MC_IP_lineages.png"))
	plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
	sapply(LineageTree[ c(MC_linId, IP_linId)], function(x) {lines(t(x)); points(t(x), pch = 16)})
dev.off()
}


