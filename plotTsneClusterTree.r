plotTsneClusterTree <- function( seuratObj, plotDir){

dir.create( file.path( plotDir, "clusterTree"), showWarnings = FALSE)
plotDir	<- file.path( plotDir, "clusterTree")

png( file.path( plotDir, paste0("SNN_Clus.png")))
	PlotClusterTree( seuratObj) 
dev.off()
}




