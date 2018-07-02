plot2Dcells <- function( PlotVals, clIdent, clustTypes, plotDir){

#This snippet plots distribution of cells ready to draw lineges. It does not call png or dev.off()  

cellColors <- rep("gray", length(plotVals[,1]))
cellColors[ which(clIdent == clustTypes[ "I"]) ]  	<- "cyan"
cellColors[ which(clIdent == clustTypes[ "M"]) ]	<- "black"
cellColors[ which(clIdent == clustTypes[ "Tl"]) ]	<- "red"	

plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
}


