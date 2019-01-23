#source("R/seuratNormImputedWT.r")

#remove mutant cells

ipmc		<- StashIdent(object = ipmc, save.name = "Initial_cellTypes")

#the folloing loop is over different comps, creates subdirectories for each number of comps

for (comps in 15){


#resDir is created by seuratNorm

#resDir <- file.path( getwd(), "Imputed")
#dir.create( resDir, showWarnings = FALSE)

compsDir 	<- file.path( seedDir, paste0( "2c", comps))
dir.create( compsDir, showWarnings = FALSE)
plotCompsDir	<- file.path( compsDir, "Plot")
dir.create( plotCompsDir, showWarnings = FALSE)
ipmc		<- RunTSNE(ipmc, genes.use = allGenes, seed.use = set.seed( randSeed))

source("R/initTypesGenesTsne")
plotInitTypesGenesTsne( ipmc, plotCompsDir) 

ipmc	 	<- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE)
ipmc	 	<- RunTSNE(ipmc, dims.use = 2:comps, theta = 0, perplexity = 15)

source("R/plotInitTypesPcaTsne.r")
plotInitTypesPcaTsne( ipmc, plotCompsDir)


for (resolDec in c(30)){ 

resolDir 	<- file.path( compsDir, paste0( "r", resolDec))
dir.create( resolDir, showWarnings = FALSE)
plotResolDir	<- file.path( resolDir, "Plot")
dir.create( plotResolDir , showWarnings = FALSE)


if(!require("e1071")){
  install.packages("e1071")
}

ipmc 		<- BuildSNN(ipmc, reduction.type = "pca", dims.use = 2:comps)
ipmc	 	<- FindClusters(ipmc, reduction.type = "pca", dims.use = 2:comps, resolution = resolDec/10, print.output = 0)
ipmc		<- StashIdent(object = ipmc, save.name = "PCA_nonValidated")

ipmc		<- ValidateClusters(ipmc, pc.use = 2:comps)
ipmc		<- StashIdent(object = ipmc, save.name = "PCA_Validated")
ipmc 		<- BuildClusterTree(ipmc, pcs.use = 2:comps, do.reorder = TRUE, reorder.numeric = TRUE, show.progress = FALSE, do.plot = FALSE)



source("R/plotTsneClusterTree.r")
plotTsneClusterTree( ipmc, plotResolDir) 
clust <- ipmc@ident

source("R/getClusterTypes.r")
clustTypes 	<- getClusterTypes(ipmc@ident)
if ( clustTypes == "Bad cluster cell types !") { 
cat( file = file.path( seedDir, "Error.txt"), "Bad cluster cell types !" ) 
next 

} 

ipmcMD 	<- RunTSNE(ipmc, dims.use = 2:comps, theta = 0, perplexity = 15, dim.embed = comps) #we need this to get lineages and principal curves with slingShot

if(!require(slingshot)){
  install.packages("slingshot")
}

tSNEValsMD	<- as.matrix(ipmcMD@dr$tsne@cell.embeddings)
slingObjMD 	<- slingshot(tSNEValsMD, ipmcMD@ident, start.clus = clustTypes["Tl"],end.clus=c(clustTypes["I"], clustTypes["M"]))

MC_linId	<- which( as.numeric( slingObjMD@slingParams$end.clus) == clustTypes["M"])
MC_linName	<- paste0("Lineage", MC_linId)
IP_linId	<- which( as.numeric( slingObjMD@slingParams$end.clus) == clustTypes["I"])
IP_linName	<- paste0("Lineage", IP_linId)

ipmc2D		<- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15, dim.embed = 2)
plotVals	<- ipmc2D@dr$tsne@cell.embeddings

source("R/getLineageCoords.r")
LineageTree	<- getLineageCoords(ipmc2D, slingObjMD)  

linPlotDir <- ( file.path( plotResolDir, "lineagePlots"))
dir.create( linPlotDir,  showWarnings = FALSE)

source("R/plot2DallLineages.r")
plot2DallLineages( LineageTree, plotVals, clustTypes, plotResolDir, IP_linId)

#the following part plot lineage plots for all lineages

source("R/lineageVlnPlot.r")
source("R/plot2DidLineage.r")
source("R/plot2Dcells.r")

#allGenes <- "sox10" # NB !!!

#	for (lineageId in seq(1, length( slingObjMD@lineages))){
#		linIdPlotDir <- file.path( linPlotDir, names( slingObjMD@lineages)[lineageId] )
#		dir.create( linIdPlotDir, showWarnings = FALSE)
#		png( file.path( linIdPlotDir, paste0( "Lineage", lineageId, "_plot.png")))
 #       		plot2Dcells( plotVals, ipmc@ident, clustTypes, linIdPlotDir)
#			plot2DidLineage( LineageTree, lineageId)
#		dev.off()	
#Now plot violin plots for all genes for the lineages
#		for (gene in allGenes) lineageVlnPlot( ipmcMD, slingObjMD, gene, lineageId, MC_linId, IP_linId, linIdPlotDir)
#	}

} #comps
} #resolDec 


cat("Now getting principle curves\n")

psTime		<- slingPseudotime(slingObjMD)
source("R/getPseudoOrder.r")
curves		<- getPseudoOrder(slingObjMD)


