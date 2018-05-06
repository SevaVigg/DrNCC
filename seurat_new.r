

#This snippet reads cells from the file Deduplicated.csv and creates _ipmc, the Seurat object. Initial cell types are found in _types variable
#new cell types _newtypes are  LETTERS A-H, W, I, M due to the preferrable single letter cell tags by Seurat visualizers.   
# Remember that the ipmc@ident stores the results of clustering, so FindCluster destroys it

require(Seurat)
require(methods)

plotDir <- file.path(getwd(), "Plot")
resDir	<- file.path(getwd(), "Res")

Genes 	<- read.table( file = paste0("Res", .Platform$file.sep, "expressionTableDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Genes_nh<- read.table( file = paste0("Res", .Platform$file.sep, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells	<- read.table( file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)
Probes	<- read.table( file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)

rownames(Genes)[which(rownames(Genes)=="Kanamycin Pos")] <- "Kanamycin_Pos"

cells_ind 	<- order(as.numeric(Cells["hpf",]))			# order with hpf increasing
Genes		<- Genes[, cells_ind]
Genes_nh	<- Genes_nh[, cells_ind]
Cells		<- Cells[, cells_ind]

types 		<- unique(paste0(Cells["hpf",], "_", Cells["CellType",]))
hpf_CellType	<- t(data.frame(hpf_CellType = paste0(Cells["hpf",], "_", Cells["CellType",]), row.names = colnames(Cells)))
Cells		<- rbind(Cells, hpf_CellType)

newTypes   	<- c("18", "21", "24", "Tl", "30", "W2", "m6", "36", "48", "I", "M", "60", "72")
names(newTypes)	<- types
dropTypes	<- c("30_standard", "36_standard", "48_standard", "60_standard", "72_standard")
types30		<- types[ which( ! types %in% dropTypes)]

Genes30		<- Genes[, which(! Cells["hpf_CellType",] %in% dropTypes)]
Genes30_nh	<- Genes_nh[, which(! Cells["hpf_CellType",] %in% dropTypes)]
Cells30		<- Cells[, which(! Cells["hpf_CellType",] %in% dropTypes)]
 

allGenes	<- rownames(Genes)
allGenes_nh	<- rownames(Genes_nh)
logExps 	<- log10(1+Genes)
logExps_nh	<- log10(1+Genes_nh)
logExps30	<- log10(1+Genes30)
logExps30_nh	<- log10(1+Genes30_nh)

ipmc    	<- CreateSeuratObject( raw.data = logExps )
ipmc_nh		<- CreateSeuratObject( raw.data = logExps_nh)
ipmc30		<- CreateSeuratObject( raw.data = logExps30)
ipmc30_nh  	<- CreateSeuratObject( raw.data = logExps30_nh)

ipmc    	<- AddMetaData( object = ipmc, t(Cells), col.name = rownames(Cells) )
ipmc_nh 	<- AddMetaData( object = ipmc_nh, t(Cells), col.name = rownames(Cells) )
ipmc30    	<- AddMetaData( object = ipmc30, t(Cells30), col.name = rownames(Cells) )
ipmc30_nh 	<- AddMetaData( object = ipmc30_nh, t(Cells30), col.name = rownames(Cells) )

newTypeDF	<- data.frame( newType = character(ncol(Cells)), row.names = colnames(Cells) )
newTypeDF30	<- data.frame( newType = character(ncol(Cells30)), row.names = colnames(Cells30) )

cellNamesDF	<- data.frame( cellNames = colnames(Cells), row.names = colnames(Cells))

ipmc		<- AddMetaData( object = ipmc, newTypeDF, col.name = "newType")
ipmc		<- AddMetaData( object = ipmc, cellNamesDF, col.name = "cellNames")
ipmc_nh		<- AddMetaData( object = ipmc_nh, newTypeDF, col.name = "newType")
ipmc_nh		<- AddMetaData( object = ipmc_nh, cellNamesDF, col.name = "cellNames")
ipmc30		<- AddMetaData( object = ipmc30, newTypeDF30, col.name = "newType")
ipmc30_nh	<- AddMetaData( object = ipmc30_nh, newTypeDF30, col.name = "newType")

levels(ipmc@ident)    	<- newTypes
levels(ipmc30@ident) 	<- newTypes[types30]
levels(ipmc_nh@ident) 	<- newTypes
levels(ipmc30_nh@ident)	<- newTypes[types30]

ipmc@ident	<- as.factor(unlist( lapply( ipmc@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
ipmc_nh@ident 	<- as.factor(unlist( lapply( ipmc_nh@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
names(ipmc@ident) <- names(ipmc_nh@ident) <- colnames(ipmc@data)

ipmc30@ident 	<- as.factor(unlist( lapply( ipmc30@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
ipmc30_nh@ident <- as.factor(unlist( lapply( ipmc30_nh@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
names(ipmc30@ident) <- names(ipmc30_nh@ident) <- colnames(ipmc30@data)

ipmc@meta.data$newType 		<- ipmc@ident
ipmc_nh@meta.data$newType	<- ipmc_nh@ident
ipmc30@meta.data$newType	<- ipmc30@ident
ipmc30_nh@meta.data$newType	<- ipmc30_nh@ident

ipmc		<- ScaleData(ipmc)
ipmc_nh		<- ScaleData(ipmc_nh)
ipmc30		<- ScaleData(ipmc30)
ipmc30_nh	<- ScaleData(ipmc30_nh)

cellColors <- paste0("gray", seq(50+3*length(levels(ipmc_nh@ident)), 50, -3))
cellColors[ which(levels(ipmc_nh@ident) == "I")] = "green3"
cellColors[ which(levels(ipmc_nh@ident) == "M")] = "black"
cellColors[ which(levels(ipmc_nh@ident) == "m6")] = "red"
cellColors[ which(levels(ipmc_nh@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc_nh@ident) == "Tl")] = "brown"

ipmc_nh <- RunTSNE(ipmc_nh, genes.use = rownames(ipmc_nh@data))
TSNEPlot(ipmc_nh, colors.use = cellColors)

mutants_keep <- which(ipmc_nh@ident %in% c("m6","W2"))



