


#This snippet identifies clusters in the PCA space and appends cluster ids to the cell annotation file
#it uses ipmc Seurat object prepared by CalcTSNE_PCASpace.r, which also supplies the number of the PCA
#components, comps

dir.create(file.path(getwd(), "Res", "PCASpaceClusters"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Plot", "PCAClusterPlot"), showWarnings = FALSE)

resDir 	<- file.path(getwd(), "Res", "PCASpaceClusters")
plotDir	<- file.path(getwd(), "Plot", "PCAClusterPlot")

#Fetch Clusters

#resolDec	<- 75			#resolution multiplied by 10, needed to mark files, if commented must be assigned
#in the calling file

ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10, print.output = 0)

png(file.path(plotDir, paste0("SNN_ClusterTreePCAS", "_c", comps, "_r", resolDec, ".png")))
	ipmc 	<- BuildClusterTree(ipmc, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE)
dev.off()





