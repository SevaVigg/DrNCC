
#This snippet identifies clusters in the gene space (using all genes) and appends cluster ids to the cell annotation file

#dir.create(file.path(getwd(), "Res"), showWarnings = FALSE)
#dir.create(file.path(getwd(), "Plot"), showWarnings = FALSE)


#dir.create(file.path(getwd(), "Plot", "GeneClusterPlot"), showWarnings = FALSE)
#plotDir	<- file.path(getwd(), "Plot", "GeneClusterPlot")


#Draw initial TSNEPlot
ipmc <- RunTSNE(ipmc, genes.use = allGenes, do.fast = FALSE)

cellColors <- paste0("gray", seq(50+3*length(levels(ipmc@ident)), 50, -3))
cellColors[ which(levels(ipmc@ident) == "I")] = "green3"
cellColors[ which(levels(ipmc@ident) == "M")] = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "red"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "brown"

png(file.path(plotDir, paste0("TSNECellsGenes", ".png")))
	TSNEPlot(ipmc, colors.use = cellColors)
dev.off()





