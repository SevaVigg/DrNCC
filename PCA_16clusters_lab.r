dir.create(file.path(getwd(), "Plot", "3dLabPlots"), showWarnings = FALSE)
plotDir	<- file.path(getwd(), "Plot", "3dLabPlots")

source("R/seuratNorm.r")

comps <- 6 

source("R/CalcTSNE_PCASpace.r")

resolDec <- 75 

source("R/FetchTSNEClustersPCASpace.r")

MIGN <- SubsetData(ipmc, ident.use = c(15:35), do.center = TRUE, do.scale = TRUE, do.clean = TRUE)
MIGN <- ScaleData(MIGN, do.scale = TRUE, do.center = TRUE)
MIGN <- RunPCA(MIGN, pc.genes = allGenes, do.print = FALSE)
MIGN <- FindClusters(MIGN, dims.use = 1:3, resolution = 0.4, prune.SNN = 0, print.output = FALSE)
MIGN <- BuildClusterTree(MIGN, pcs.use = 1:3, do.reorder = TRUE, reorder.numeric = TRUE)

cellColors[1]<- "black"
cellColors[2]<- "red"
cellColors[3]<- "green"
cellColors[4]<-  "grey"
cellColors[5] <- "orange"

png(file.path(plotDir, paste0("PCAPlot3d_12", ".png")))
	PCAPlot(MIGN, dim.1 = 1, dim.2 = 2, cols.use = cellColors)
dev.off()

png(file.path(plotDir, paste0("PCAPlot3d_13", ".png")))
	PCAPlot(MIGN, dim.1 = 1, dim.2 = 3, cols.use = cellColors)
dev.off()

png(file.path(plotDir, paste0("PCAPlot3d_23", ".png")))
	PCAPlot(MIGN, dim.1 = 2, dim.2 = 3, cols.use = cellColors)
dev.off()

MIGN <- RunTSNE(MIGN, dims.use = 1:3, do.fast = FALSE, dim.embed = 3, perplexity = 100)

png(file.path(plotDir, paste0("TSNEPlot3d_12", ".png")))
	TSNEPlot(MIGN, dim.1 = 1, dim.2 = 2, colors.use = cellColors)
dev.off()

png(file.path(plotDir, paste0("TSNEPlot3d_13", ".png")))
	TSNEPlot(MIGN, dim.1 = 1, dim.2 = 3, colors.use = cellColors)
dev.off()

png(file.path(plotDir, paste0("TSNEPlot3d_23", ".png")))
	TSNEPlot(MIGN, dim.1 = 2, dim.2 = 3, colors.use = cellColors)
dev.off()




