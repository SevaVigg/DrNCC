#This snippet calculates tSNE of ipmc object (using comps PCA components) 

dir.create(file.path(getwd(), "Plot", "PCAClusterPlot"), showWarnings = FALSE)
plotDir	<- file.path(getwd(), "Plot", "PCAClusterPlot")

ipmc <- RunPCA(object = ipmc, pc.genes = allGenes, pcs.compute = 20, do.print = FALSE)

#Draw initial TSNEPlot
ipmc <- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15)

initCellColors <- rep( "grey", length(levels(ipmc@ident)))
initCellColors[ levels(ipmc@ident) == "I"] <- "cyan"
initCellColors[ levels(ipmc@ident) == "M"] <- "black"
initCellColors[ levels(ipmc@ident) == "Tl"] <- "red"
initCellColors[ levels(ipmc@ident) == "m6"] <- "blue"
initCellColors[ levels(ipmc@ident) == "W2"] <- "magenta"

png(file.path(plotDir, paste0("TSNE_CellsPCAS", "_c", comps, ".png")))
	TSNEPlot(ipmc, colors.use = initCellColors)
dev.off()


