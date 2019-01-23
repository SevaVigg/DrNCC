plotInitTypesGenesTsne <- function( seuratObj, plotDir){

dir.create( file.path( plotDir, "InitTypes_Genes_TSNE"), showWarnings = FALSE)
plotDir	<- file.path( plotDir, "InitTypes_Genes_TSNE")

#seuratObj  <- RunTSNE( seuratObj, genes.use = allGenes, seed.use = as.numeric(as.POSIXct(Sys.time())), tsne.method = "tsne", perplexity = 20)
seuratObj  <- RunTSNE( seuratObj, genes.use = allGenes, seed.use = as.numeric(as.POSIXct(Sys.time())), theta = 0, 
			eta = 10, max_iter = 5000, 
			perplexity = 20, verbose = FALSE)

cellColors <- rev( sequential_hcl( length( levels(seuratObj@ident)), h = 100, c. = c(60, 100), l = c(60, 100)))

cellColors[ which(levels(ipmc@ident) == "I")]  = "cyan"
cellColors[ which(levels(ipmc@ident) == "M")]  = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "pink"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "red"

png( file.path( plotDir, "initCellsTsne.png"))
	TSNEPlot(seuratObj, colors.use = cellColors)
dev.off()
}

