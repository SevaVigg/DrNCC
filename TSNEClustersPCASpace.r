
#This snippet identifies clusters in the PCA space and appends cluster ids to the cell annotation file
#it uses ipmc Seurat object prepared by CalcTSNE_PCASpace.r, which also supplies the number of the PCA
#components, comps

#dir.create(file.path(getwd(), "Res"), showWarnings = FALSE)
#dir.create(file.path(getwd(), "Plot"), showWarnings = FALSE)

dir.create(file.path(getwd(), "Res", "PCASpaceClusters"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Plot", "PCAClusterPlot"), showWarnings = FALSE)

resDir 	<- file.path(getwd(), "Res", "PCASpaceClusters")
plotDir	<- file.path(getwd(), "Plot", "PCAClusterPlot")

#Fetch Clusters

resolDec	<- 6			#resolution multiplied by 10, needed to mark files
ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10, print.output = 0)

png(file.path(plotDir, paste0("SNN_ClusterTreePCAS", "_c", comps, "_r", resolDec, ".png")))
	ipmc 	<- BuildClusterTree(ipmc, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE)
dev.off()

Cells <- rbind(Cells, ipmc@ident)
rownames(Cells)[length(rownames(Cells))]	<- paste0("ClusterPCASp", "_c", comps, "_r", resolDec)

png(file.path(plotDir, paste0("TSNE_ClustersPCAS", "_c", comps, "_r", resolDec, ".png")))
	TSNEPlot(ipmc)
dev.off()

#write.table(Cells, file = file.path( resDir, "CellClustersAll.tsv"), sep = "\t")
#write.table(Genes, file = file.path( resDir, "GenesAllCl.tsv"), sep = "\t")







