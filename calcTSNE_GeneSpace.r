#Draw initial TSNEPlot

calcTSNE_GeneSpace <- function( seuratObj){

source("R/setCellTypeColors.r")

seuratObj <- RunTSNE( seuratObj, genes.use = rownames(seuratObj@data), seed.use = as.numeric(as.POSIXct(Sys.time())), theta = 0, 
	eta = 10, max_iter = 5000, 
	perplexity = 30, verbose = FALSE)

	TSNEPlot(seuratObj, colors.use = setCellTypeColors(seuratObj))

return( seuratObj)
}




