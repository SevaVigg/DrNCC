#requires logExp matrix

perplexityTestDir <- file.path( resDir, "perplexityTest")
dir.create( perplexityTestDir, showWarnings = FALSE)

for (seed in 31:33){

seedDir	<- file.path( perplexityTestDir, paste0("s", seed))
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

source("R/setCellTypeColors.r")
for (perplexity in c(2, 5, 10, 15, 20, 30, 40, 50, 75, 100, 150)){
	cat("perplexity = ", perplexity, "\n")
	perpDir 	<- file.path( seedDir, paste0("p", perplexity))
	dir.create( perpDir, showWarnings = FALSE)
		for (iter in 1:5){
		seedFile	<- file.path( perpDir, paste0("seed_", iter))
		seedPerp	<- as.numeric(as.POSIXct(Sys.time()))		 	
		cat(file = seedFile, seedPerp, "\n")
		ipmcPerp  <- RunTSNE( ipmc, genes.use = allGenes, seed.use = seedPerp, theta = 0, 
			eta = 10, max_iter = 3000, 
			perplexity = perplexity, verbose = FALSE)
		perpFile 	<- file.path( perpDir, paste0("p", perplexity, "_", iter, ".png"))
		png( perpFile) 
			TSNEPlot( ipmcPerp, colors.use = setCellTypeColors(ipmcPerp))
		dev.off()
		rm(ipmcPerp)
		close(seedFile)
		}	#iteration
	} #perplexity
} #seed		

