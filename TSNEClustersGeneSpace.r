#This snippet identifies clusters in the gene space (using all genes) and appends cluster ids to the cell annotation file

#dir.create(file.path(getwd(), "Res"), showWarnings = FALSE)
#dir.create(file.path(getwd(), "Plot"), showWarnings = FALSE)

dir.create(file.path(getwd(), "Res", "geneSpaceClusters"), showWarnings = FALSE)
dir.create(file.path(getwd(), "Plot", "ClusterPlot"), showWarnings = FALSE)

resDir 	<- file.path(getwd(), "Res", "geneSpaceClusters")
plotDir	<- file.path(getwd(), "Plot", "ClusterPlot")


#Draw initial TSNEPlot
cellColors <- paste0("gray", seq(50+3*length(levels(ipmc@ident)), 50, -3))
cellColors[ which(levels(ipmc@ident) == "I")] = "green3"
cellColors[ which(levels(ipmc@ident) == "M")] = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "red"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "brown"

ipmc <- RunTSNE(ipmc, genes.use = rownames(ipmc@data), do.fast = FALSE)
TSNEPlot(ipmc, colors.use = cellColors)

#Fetch Clusters
params		<- 3 
ipmc 	<- FindClusters(ipmc, genes.use = rownames(ipmc@data), k.param = params, k.scale = 1000)

png(file.path(plotDir, paste0("SNNClusterTreeGSallGenes", "_k", params, ".png")))
	ipmc 	<- BuildClusterTree(ipmc, genes.use = allGenes, do.reorder = TRUE, reorder.numeric = TRUE)
dev.off()

#Cells <- rbind(Cells, ipmc@ident)
#rownames(Cells)[length(rownames(Cells))]	<- "Cluster_geneSp"

png(file.path(plotDir, paste0("TSNEClustersGSallGenes", "_k", params, ".png")))
	TSNEPlot(ipmc)
dev.off()

write.table(Cells, file = file.path( resDir, "CellClustersAll.tsv"), sep = "\t")
write.table(Genes, file = file.path( resDir, "GenesAllCl.tsv"), sep = "\t")







