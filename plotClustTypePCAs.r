plotClustTypePCAs <- function( seuratObj, clustTypes, plotDir, Ncomps){

#This snippet plots distribution of cells ready to draw lineges. It does not call png or dev.off()  

PCAPlotDirName 	<- file.path(plotDir, "PCAPlots")
dir.create(PCAPlotDirName, showWarnings = FALSE)

clIdent <- seuratObj@ident

cellColors <- rep("gray", length( unique(clIdent)))
cellColors[ clustTypes[ "I"] ]  	<- "cyan"
cellColors[ clustTypes[ "M"] ]	<- "black"
cellColors[ clustTypes[ "Tl"] ]	<- "red"	

for (c1 in 1:(Ncomps-1)){
	for (c2 in (c1+1):Ncomps){
		png( file.path( PCAPlotDirName, paste0( "PCA_", c1, "_", c2, ".png")))
		PCAPlot( seuratObj, dim.1 = c1, dim.2 = c2, cols.use = cellColors)
		dev.off()
	}

}
   
}


