# This script identifies trajectories through clusters obtained by 
# cell subsampling 

source("R/seuratNorm.r")
resDir		<-	file.path( getwd(), "Res")
subSDir 	<- 	file.path( resDir, "Subsampling")
dir.create( subSDir, showWarnings = FALSE)
    
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
ipmc2D		<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = 2)
plotVals	<- ipmc2D@dr$tsne@cell.embeddings

cat("now calculate MD tSNE", "\n")
#ipmcMD  	<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = comps) 
#we need this to get lineages and principal curves with slingShot
#tSNEValsMD      <- as.matrix(ipmcMD@dr$tsne@cell.embeddings)

source("R/plotInitTypesPcaTsne.r")
#plotInitTypesPcaTsne( ipmc2D, compsDir)

resolDec 	<- 80

resolDir 	<- file.path( compsDir, paste0( "r", resolDec))
dir.create( resolDir, showWarnings = FALSE)
plotResolDir	<- file.path( resolDir, "Plot")
dir.create( plotResolDir , showWarnings = FALSE)

ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10, print.output = 0)
ipmc 		<- BuildClusterTree(ipmc, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE, do.plot = FALSE)

source("R/plotTsneClusterTree.r")
plotTsneClusterTree( ipmc, plotResolDir) 


#NO BRANCH
#removing the branch
#	noBranch 	<- WhichCells(ipmc, ident = c(29:53))
#	ipmcNB 		<- SubsetData(ipmc, ident.remove = c(29:53))
#	noBranchDir	<- file.path( resolDir, "noBranch")
#	dir.create( noBranchDir, showWarnings = FALSE)
#	plotResolDir	<- noBranchDir
#
#	# calculate PCA
#	cat("now calculate noBranch PCA", "\n") 
#	ipmcNB          <- RunPCA(ipmcNB, pc.genes = rownames( ipmcNB@data), pcs.compute = comps, do.print = FALSE) 
#	# this is 2D tSNE for all cells, which we will use for drawing
#	cat("now calculate 2D tSNE", "\n")
#	ipmc2DNB	<- RunTSNE(ipmcNB, dims.use = 1:comps, theta = 0, perplexity = 15)
#	plotVals	<- ipmc2DNB@dr$tsne@cell.embeddings
#
#BRANCH REMOVAL PROCEDURE ENDs

source("R/getClusterTypes.r")
allClusterTypes <- getClusterTypes( ipmc@ident)

if(!require(slingshot)){
  install.packages("slingshot")
}
source("R/getLineageCoords.r")
source("R/plot2DallLineages.r")
source("R/plot2DidLineage.r")
source("R/plot2Dcells.r")

png( file.path( plotResolDir, "IP_subsampling.png"))
	plot2Dcells( plotVals, ipmc@ident, allClusterTypes, plotResolDir)  #Note that only ipmc contains correct clustering, not ipmc2D

MC_linList	<- list()
MC_linTree	<- numeric(0)
IP_linList	<- list()
IP_linTree	<- numeric(0)

MC_linFileName	<- file.path( resolDir, "MC_lineages.txt")
IP_linFileName	<- file.path( resolDir, "IP_lineages.txt")
sampleFileName	<- file.path( resolDir, "sample.txt")
MC_linFile 	<- file( MC_linFileName, open = "w")
IP_linFile	<- file( IP_linFileName, open = "w")
sampleFile	<- file( sampleFileName, open = "w")

nRounds		<- 5
keepShare	<- 0.9


sampleMatrix	<- sapply(1:nRounds, function(x) sample(ipmc@cell.names, size = round( length( ipmc@cell.names)*keepShare), replace = FALSE))
write.table(sampleMatrix, file = sampleFile)
close( sampleFile)

for ( subRound in 1:nRounds ){
	subSample	<- sampleMatrix[ , subRound ]
	ipmcSub		<- SubsetData( ipmc, cells.use = subSample) #this is a SeuratObject containing only cells in the subRound
	ipmcSub		<- RunPCA( ipmcSub, pc.genes = allGenes, pcs.compute = comps, do.print = FALSE) 
	ipmcSub		<- FindClusters( ipmcSub, reduction.type = "pca", dims.use = 1:comps, resolution = resolDec/10,  print.output = 0, force.recalc = TRUE) 
	ipmcSub 	<- BuildClusterTree(ipmcSub, genes.use = allGenes, pcs.use = 1:comps, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = FALSE, show.progress = FALSE)
	ipmcSubMD  	<- RunTSNE( ipmcSub, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = comps) 
	#we need this to get lineages and principal curves with slingShot
	tSNEValsMD      <- as.matrix( ipmcSubMD@dr$tsne@cell.embeddings)
	clustTypes	<- getClusterTypes( ipmcSub@ident)	
	slingObjRound 	<- slingshot( tSNEValsMD, ipmcSub@ident, start.clus = clustTypes["Tl"], end.clus = c(clustTypes["I"], clustTypes["M"]) )
	MC_linId	<- which( as.numeric( slingObjRound@lineageControl$end.clus) == clustTypes["M"])
	MC_linList[[subRound]]	<- slingObjRound@lineages[[MC_linId]] 
	MC_linName	<- paste0("Lineage", MC_linId)
	IP_linId	<- which( as.numeric( slingObjRound@lineageControl$end.clus) == clustTypes["I"])
	IP_linName	<- paste0("Lineage", IP_linId)
	IP_linList[[subRound]]	<- slingObjRound@lineages[[IP_linId]]
	LineageTree	<- getLineageCoords( ipmc2D, slingObjRound) 
	lineageId	<- MC_linId 
	plot2DidLineage( LineageTree, lineageId)	
	lineageId	<- IP_linId 
	plot2DidLineage( LineageTree, lineageId)	
}
dev.off()
lapply(MC_linList, function(x) cat( file = MC_linFile, x, "\n"))
lapply(IP_linList, function(x) cat( file = IP_linFile, x, "\n"))
close( MC_linFile)
close( IP_linFile)
