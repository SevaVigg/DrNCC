plotAllLineagePlots <- function( slingObj, LineageTree, plotDir ){
	source("R/plot2Dcells.r")
	source("R/plot2DidLineage.r")
	for (lineageId in seq(1, length( slingObj@lineages))){
		linIdPlotDir <- file.path( plotDir, names( slingObj@lineages)[lineageId] )
		dir.create( linIdPlotDir, showWarnings = FALSE)
		png( file.path( linIdPlotDir, paste0( "Lineage", lineageId, "_plot.png")))
        		plot2Dcells( plotVals, ipmc@ident, clustTypes, linIdPlotDir)
			plot2DidLineage( LineageTree, lineageId)
		dev.off()	
		for (gene in allGenes) lineageVlnPlot( ipmcMD, slingObj, gene, lineageId, MC_linId, IP_linId, linIdPlotDir)
	}
