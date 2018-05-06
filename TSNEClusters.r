resDir <- paste0(getwd(), .Platform$file.sep, "Res")

ipmc_orig	<- ipmc_nh

ipmc_nh 	<- FindClusters(ipmc_nh, genes.use = rownames(ipmc_nh@data), k.param = 3, k.scale = 1000)

png("Plot/SNNClusterTree.png")
	ipmc_nh 	<- BuildClusterTree(ipmc_nh, genes.use = allGenes, do.reorder = TRUE, reorder.numeric = TRUE)
dev.off()

colnames(Genes_nh) <- paste0( ipmc_nh@ident, "_", colnames(Genes_nh))
colnames(Cells)	   <- paste0( ipmc_nh@ident, "_", colnames(Cells))

ipmc_nh <- RunTSNE(ipmc_nh, genes.use = allGenes_nh)

Cells <- rbind(Cells, ipmc_nh@ident)
rownames(Cells)[length(rownames(Cells))]	<- "Cluster_ids"

png("Plot/TSNEClustersAll.png")
TSNEPlot(ipmc_nh)
dev.off()

write.table(Cells, file = file.path( resDir, "CellClustersAll.tsv"), sep = "\t")
write.table(Genes_nh, file = file.path( resDir, "GenesAllCl.tsv"), sep = "\t")





