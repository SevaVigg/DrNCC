source("R/seuratNorm.r")

origIdent 	<- ipmc@ident

ipmc		<- SubsetData(ipmc, ident.remove = c("W2", "m6"))

#for (comps in seq(5, 20, 5)){


#resDir is created by seuratNorm
dir.create( file.path( resDir, paste0( "c", comps)), showWarnings = FALSE)
compsDir 	<- file.path( resDir, paste0( "c", comps))
dir.create( file.path( compsDir, "Plot" ) , showWarnings = FALSE)
plotCompsDir	<- file.path( compsDir, "Plot")

ipmc	 	<- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE)
ipmc	 	<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15)
source("R/plotInitTypesPcaTsne.r")
plotInitTypesPcaTsne( ipmc, plotCompsDir)


for (resolDec in seq(20, 100, 20)){ 

dir.create( file.path( compsDir, paste0( "r", resolDec)), showWarnings = FALSE)
resolDir 	<- file.path( compsDir, paste0( "r", resolDec))
dir.create( file.path( resolDir, "Plot" ) , showWarnings = FALSE)
plotResolDir	<- file.path( resolDir, "Plot")

ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10, print.output = 0)
ipmc 		<- BuildClusterTree(ipmc, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE, do.plot = FALSE)
source("R/plotTsneClusterTree.r")
plotTsneClusterTree( ipmc, plotResolDir) 
clust <- ipmc@ident

source("R/getClusterTypes.r")
clustTypes 	<- getClusterTypes(ipmc@ident)

ipmcMD 	<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = comps) #we need this to get lineages and principal curves with slingShot

if(!require(slingshot)){
  install.packages("slingshot")
}

tSNEValsMD	<- as.matrix(ipmcMD@dr$tsne@cell.embeddings)
slingObjMD 	<- slingshot(tSNEValsMD, clust, start.clus = clustTypes["Tl"],end.clus=c(clustTypes["I"], clustTypes["M"]))

MC_linId	<- which( as.numeric( slingObjMD@lineageControl$end.clus) == clustTypes["M"])
MC_linName	<- paste0("Lineage", MC_linId)
IP_linId	<- which( as.numeric( slingObjMD@lineageControl$end.clus) == clustTypes["I"])
IP_linName	<- paste0("Lineage", IP_linId)

ipmc_2D		<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = 2)
plotVals	<- ipmc_2D@dr$tsne@cell.embeddings

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords(ipmc_2D, slingObjMD)  

dir.create( file.path( plotResolDir, "lineagePlots"), showWarnings = FALSE)
linPlotDir <- ( file.path( plotResolDir, "lineagePlots"))

source("R/plot2DallLineages.r")
plot2DallLineages( LineageTree, plotVals, clustTypes, plotResolDir)

source("R/lineageVlnPlot.r")
source("R/plot2DidLineage.r")

#allGenes <- "sox10" # NB !!!

for (lineageId in seq(1, length( slingObjMD@lineages))){
	dir.create( file.path( linPlotDir, names(slingObjMD@lineages)[lineageId] ), showWarnings = FALSE)
	linIdPlotDir <- file.path( linPlotDir, names(slingObjMD@lineages)[lineageId] )
	plot2DidLineage(ipmcMD, LineageTree, lineageId, MC_linId, IP_linId, plotVals, clustTypes, linIdPlotDir)	
	for (gene in allGenes) lineageVlnPlot( ipmcMD, slingObjMD, gene, lineageId, MC_linId, IP_linId, linIdPlotDir)
}
} #comps
} #resolDec 


cat("Now getting principle curves\n")

slingObjMD	<- getCurves(slingObjMD, extend = "n", shrink = TRUE, thresh = 0.001, maxit = 10)
psTime		<- pseudotime(slingObjMD)
source("R/getPseudoOrder.r")
curves		<- getPseudoOrder(slingObjMD)


