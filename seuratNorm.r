
#This is a minimal version of seurat preparing

#This snippet reads cells from the file Deduplicated.csv and creates _ipmc, the Seurat object. Initial cell types are found in _types variable
#new cell types _newtypes are  LETTERS A-H, W, I, M due to the preferrable single letter cell tags by Seurat visualizers.   
# Remember that the ipmc@ident stores the results of clustering, so FindCluster destroys it

require(Seurat)
require(methods)

#plotDir <- file.path(getwd(), "Plot")
resDir	<- file.path(getwd(), "Res")

Genes	<- read.table( file = paste0("Res", .Platform$file.sep, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells	<- read.table( file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)
Probes	<- read.table( file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)

rownames(Genes)[which(rownames(Genes)=="Kanamycin Pos")] <- "Kanamycin_Pos"

#reorder cells
cells_ind 	<- order(as.numeric(Cells["hpf",]))			# order with hpf increasing
Genes		<- Genes[, cells_ind]
Cells		<- Cells[, cells_ind]

#rename cell types, prepare the annotated cell table
types 		<- unique(paste0(Cells["hpf",], "_", Cells["CellType",]))
hpf_CellType	<- t(data.frame(hpf_CellType = paste0(Cells["hpf",], "_", Cells["CellType",]), row.names = colnames(Cells)))
Cells		<- rbind(Cells, hpf_CellType)

newTypes   	<- c("18", "21", "24", "Tl", "30", "W2", "m6", "36", "48", "I", "M", "60", "72")
names(newTypes)	<- types

allGenes	<- rownames(Genes)
logExps 	<- log2(1+Genes)

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExps) )
ipmc    	<- AddMetaData( object = ipmc, t(Cells), col.name = rownames(Cells) )

newTypeDF	<- data.frame( newType = character(ncol(Cells)), row.names = colnames(Cells) )
cellNamesDF	<- data.frame( cellNames = colnames(Cells), row.names = colnames(Cells))

ipmc		<- AddMetaData( object = ipmc, newTypeDF, col.name = "newType")
ipmc		<- AddMetaData( object = ipmc, cellNamesDF, col.name = "cellNames")

levels(ipmc@ident)    	<- newTypes

ipmc@ident		<- as.factor(unlist( lapply( ipmc@meta.data[ , "hpf_CellType"], function(cell) newTypes[as.character(cell)]) ))
names(ipmc@ident) 	<- names(ipmc@ident) <- colnames(ipmc@data)
ipmc@meta.data$newType 	<- ipmc@ident

ipmc			<- ScaleData(ipmc, do.scale = TRUE, do.center = TRUE)




