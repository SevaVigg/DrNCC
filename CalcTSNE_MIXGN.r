
#This snippet calculates tSNE (using comps PCA components)

dir.create(file.path(getwd(), "Plot", "PCAClusterPlot"), showWarnings = FALSE)
plotDir	<- file.path(getwd(), "Plot", "PCAClusterPlot")

ipmc <- RunPCA(object = ipmc, pc.genes = allGenes, pcs.compute = 20, do.print = FALSE)

#Draw initial TSNEPlot
comps		<- 5 

ipmc <- RunTSNE(ipmc, dims.use = 1:comps, do.fast = FALSE)

cellColors <- paste0("gray", seq(50+3*length(levels(ipmc@ident)), 50, -3))
cellColors[ which(levels(ipmc@ident) == "I")] = "green3"
cellColors[ which(levels(ipmc@ident) == "M")] = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "red"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "brown"

png(file.path(plotDir, paste0("TSNE_CellsPCAS", "_c", comps, ".png")))
	TSNEPlot(ipmc, colors.use = cellColors)
dev.off()


