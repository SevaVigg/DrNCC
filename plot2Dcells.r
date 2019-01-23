plot2Dcells <- function( PlotVals, cellColors){

#This snippet plots distribution of cells ready to draw lineges. It does not call png or dev.off()  
	
plot( plotVals[,2]~plotVals[,1], cex = 0.7, pch = 16, col = cellColors, xlab = "tSNE1", ylab = "tSNE2")
}


