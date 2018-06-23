plotInitTypesPcaTsne <- function( seuratObj, plotDir){

dir.create( file.path( plotDir, "InitTypes_PCA_TSNE"), showWarnings = FALSE)
plotDir	<- file.path( plotDir, "InitTypes_PCA_TSNE")

cellColors <- gray.colors( length( levels(seuratObj@ident)))

cellColors[ which(levels(ipmc@ident) == "I")]  = "cyan"
cellColors[ which(levels(ipmc@ident) == "M")]  = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "pink"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "red"

png( file.path( plotDir, "initCellsTsne.png"))
	TSNEPlot(seuratObj, colors.use = cellColors)
dev.off()
}

