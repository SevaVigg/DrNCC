# This script identifies trajectories through clusters obtained by 
# cell subsampling 

source("R/seuratNorm.r")
subSDir 	<- 	file.path( resDir, "Subsampling")
dir.create( subSDir, showWarnings = FALSE)

sourceFileName 	<- file.path( getwd(), "Subsampling", "15_samples_4.txt")
subS		<- read.table( sourceFileName, row.names = 1, stringsAsFactors = TRUE )
nClust		<- apply(subS, 2, max)
cat(nClust, "\n")

    
origIdent       <- ipmc@ident

ipmc            <- SubsetData(ipmc, ident.remove = c("W2", "m6"))

comps		<- 	15

compsDir 	<- file.path( subSDir, paste0( "c", comps))
dir.create( compsDir, showWarnings = FALSE)
plotCompsDir	<- file.path( compsDir, "Plot")
dir.create( plotCompsDir, showWarnings = FALSE)

# calculate PCA
cat("now calculate PCA", "\n") 
ipmc            <- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE) 
# this is 2D tSNE for all cells, which we will use for drawing
cat("now calculate 2D tSNE", "\n")
ipmc2D		<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15)
plotVals	<- ipmc2D@dr$tsne@cell.embeddings

#and this is MD tSNE for all cells, which we will use for embedding trajectories
cat("now calculate MD tSNE", "\n")
ipmcMD  	<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = comps) #we need this to get lineages and principal curves with slingShot
tSNEValsMD      <- as.matrix(ipmcMD@dr$tsne@cell.embeddings)

source("R/plotInitTypesPcaTsne.r")
plotInitTypesPcaTsne( ipmc2D, subSDir)

resolDec 	<- 80

resolDir 	<- file.path( compsDir, paste0( "r", resolDec))
dir.create( resolDir, showWarnings = FALSE)
plotResolDir	<- file.path( resolDir, "Plot")
dir.create( plotResolDir , showWarnings = FALSE)

ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10, print.output = 0)
ipmc 		<- BuildClusterTree(ipmc, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE, do.plot = FALSE)
source("R/plotTsneClusterTree.r")
plotTsneClusterTree( ipmc, plotResolDir) 


source("R/getClusterTypes.r")
allClusterTypes <- getClusterTypes( ipmc@ident)

if(!require(slingshot)){
  install.packages("slingshot")
}
source("R/getLineageCoords.r")
source("R/plot2DallLineages.r")
source("R/plot2DidLineage.r")
source("R/plot2Dcells.r")

png( file.path( plotResolDir, "Lineage_plot.png"))
	plot2Dcells( plotVals, ipmc@ident, allClusterTypes, plotResolDir) 

MC_linMatrix	<- integer(0)
MC_linTree	<- numeric(0)
IP_linMatrix	<- integer(0)
IP_linTree	<- numeric(0)

	#for ( subRound in 1:1 ){
	for ( subRound in 1:ncol(subS)){
	clFactor	<- factor( subS[, subRound])
	names(clFactor)	<- rownames(subS)
	clustTypes	<- getClusterTypes( clFactor )
	cellsRoundObj	<- SubsetData(ipmc, cells.use = rownames(subS)) #this is a SeuratObject containing only cells in the subRound
	tSNEValsRound	<- tSNEValsMD[rownames(subS), ]	#this is a stub; to be substituted for the round specific cells
	slingObjRound 	<- slingshot( tSNEValsRound, clFactor, start.clus = clustTypes["Tl"], end.clus = c(clustTypes["I"], clustTypes["M"]) )
	MC_linId	<- which( as.numeric( slingObjRound@lineageControl$end.clus) == clustTypes["M"])
	MC_linMatrix	<- cbind( MC_linMatrix, slingObjRound@lineages[[MC_linId]]) 
	MC_linName	<- paste0("Lineage", MC_linId)
	IP_linId	<- which( as.numeric( slingObjRound@lineageControl$end.clus) == clustTypes["I"])
	IP_linName	<- paste0("Lineage", IP_linId)
	IP_linMatrix	<- cbind(IP_linMatrix, slingObjRound@lineages[[IP_linID]])
	LineageTree	<- getLineageCoords( ipmc2D, slingObjRound) 
	#MC_linTree	<- cbind( MC_linTree, LineageTree[
	lineageId	<- MC_linId 
	plot2DidLineage( LineageTree, lineageId)	
	lineageId	<- IP_linId 
	plot2DidLineage( LineageTree, lineageId)	
}
dev.off()




clustTypes      <- getClusterTypes(ipmc@ident)
