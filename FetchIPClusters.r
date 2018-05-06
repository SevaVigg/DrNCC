# This snippet runs SNN clustering tool for _ipmc, the Seurat object, using IPMarkers prepared by getIPmarkers.r

require("Seurat")

ipmc <- FindClusters(ipmc, genes.use = names(IPmarkers), k.param = 4, k.scale = 1000, prune.SNN = 1/15, resolution = 0.1, print.output = F)

ClRes <- as.data.table(ipmc@data.info)

MaxIP_cl <- ClRes[,.(N_IP = sum(celltype == "I")),res.0.1][order(-N_IP)][1, res.0.1]

IPclGenes <- ipmc@data[,which(ipmc@data.info$"res.0.1" == MaxIP_cl)]

OutFileName 	<- file.path(getwd(),"Source","Cells","IP_friendly.txt")
OutFile 	<- file(OutFileName, open = "w")

cat(names(IPclGenes),file=OutFile)

close(OutFile)


