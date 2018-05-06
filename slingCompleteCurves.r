
require(slingshot)

seurObj <- ipmc

tSNEVals<-seurObj@dr$tsne@cell.embeddings
clust <- seurObj@ident

CellClusters <- sapply(unique(levels(ipmc@ident)), function(x) table(unlist(Cells["hpf_CellType" , WhichCells(ipmc, x)])))

cl_IP	<- which.max(sapply(CellClusters, function(x) if(length(x[grep("IP", names(x))])==0) 0 else x[grep("IP", names(x))]))
cl_MC	<- which.max(sapply(CellClusters, function(x) if(length(x[grep("MC", names(x))])==0) 0 else x[grep("MC", names(x))]))
cl_tail	<- which.max(sapply(CellClusters, function(x) if(length(x[grep("tails", names(x))])==0) 0 else x[grep("tails", names(x))]))

sds <- slingshot(tSNEVals, clust, start.clus = cl_tail ,end.clus=c( cl_IP, cl_MC ),extend="n",shrink=FALSE,thresh=0.1)

cellColors <- rep("grey", length(clust))

cellColors[as.numeric(levels(clust))[clust] == cl_MC] 	<- "black"
cellColors[as.numeric(levels(clust))[clust] == cl_tail] 	<- "red"
cellColors[as.numeric(levels(clust))[clust] == cl_IP] 	<- "cyan"
cellColors[as.numeric(levels(clust))[clust] == 17]		<- "pink"
cellColors[as.numeric(levels(clust))[clust] == 33]		<- "blue"
cellColors[as.numeric(levels(clust))[clust] == 18]		<- "magenta"

tSNE1	<- tSNEVals[ , 1]
tSNE2	<- tSNEVals[ , 2]
#tSNE3	<- tSNEVals[ , 3] 

plot(tSNE2 ~ tSNE1, col = cellColors, cex = 0.7, pch = 16)
lines(sds, lwd = 1, type="lineages", cex = 0.5)

#png(find.path(plotDir, "principal_curves.version_4.png"))
#	plot(tSNE3 ~ tSNE1, col = cellColors, cex = 0.7, pch = 16)
#	lines(sds, lwd = 1, type="curves")
#dev.off()

#png(find.path(plotDir, "lineages_version_4.png"))
#	plot(tSNE1 ~ tSNE3, col = cellColors, cex = 0.7, pch = 16)
#	lines(sds, lwd = 2,type="lineages")
#dev.off()
