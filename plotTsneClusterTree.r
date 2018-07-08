plotTsneClusterTree <- function( seuratObj, plotDir){

dir.create( file.path( plotDir, "clusterTree"), showWarnings = FALSE)
plotDir	<- file.path( plotDir, "clusterTree")

png( file.path( plotDir, paste0("SNN_Clus.png")), width = 1280, height  = 800, units = "px")
	PlotClusterTree( seuratObj) 
dev.off()
}




