
#requires logExp matrix

require(colorspace)

perplexityPCATestDir <- file.path( resDir, "perplexityPCATest")
dir.create( perplexityPCATestDir, showWarnings = FALSE)


cellColors <- rev( sequential_hcl( length( levels( ipmc@ident)), h = 100, c. = c(60, 100), l = c(60, 100)))

cellColors[ which(levels(ipmc@ident) == "I")]  = "cyan"
cellColors[ which(levels(ipmc@ident) == "M")]  = "black"
cellColors[ which(levels(ipmc@ident) == "m6")] = "pink"
cellColors[ which(levels(ipmc@ident) == "W2")] = "magenta"
cellColors[ which(levels(ipmc@ident) == "Tl")] = "red"


for (seed in 1:10){

seedDir	<- file.path( perplexityPCATestDir, paste0("s", seed))
dir.create( seedDir, showWarnings = FALSE)

randSeed <- as.integer(Sys.time())
set.seed( randSeed )
cat( file = file.path(seedDir, "seed.txt"), randSeed)
cat("seed = ", seed, "\n")

logExpsImp <- DrImpute(logExps)

#seurat scales the data by mean and sd. Let us scale the data with 

ipmc    	<- CreateSeuratObject( raw.data = as.matrix(logExpsImp) )
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
ipmc			<- StashIdent(object = ipmc, save.name = "originalCellTypes")

cat( "\n")

for (comps in c(5, 6, 7, 8, 10, 12, 15)){

compsDir 	<- file.path( seedDir, paste0( "c", comps))
dir.create( compsDir, showWarnings = FALSE)
cat("comps = ", comps, "\n")

ipmcPCA	 	<- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE)


for (perplexity in c(2, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150)){
	cat("perplexity = ", perplexity, "\n")
	perpDir 	<- file.path( compsDir, paste0("p", perplexity))
	dir.create( perpDir, showWarnings = FALSE)
		for (iter in 1:5){

		ipmcPerp  <- RunTSNE(ipmcPCA, dims.use = 1:comps, seed.use = as.numeric(as.POSIXct(Sys.time())),  theta = 0, 
			eta = 10, max_iter = 3000, 
			perplexity = perplexity, verbose = FALSE)
		perpFile 	<- file.path( perpDir, paste0("p", perplexity, "_", iter, ".png"))
		png( perpFile) 
			TSNEPlot( ipmcPerp, colors.use = cellColors)
		dev.off()
		rm(ipmcPerp)
		}	#iteration

	} #perplexity
rm(ipmcPCA)
	}#comps	
} #seed		

