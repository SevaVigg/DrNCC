lineageVlnPlot <- function(seuratObj, slingObj, gene, lineageId, MC_linId, IP_linId, plotDir){
	
dir.create( file.path( plotDir, "Vln"), showWarnings = FALSE)

source("R/setLineageColors.r")

lineage 	<- slingObj@lineages[[lineageId]]
lind 		<- data.frame(l = lineage, sl = sort( as.numeric (lineage)), stringsAsFactors = FALSE) # create a transposition
colorInd 	<- sapply( lind$l, function(x) which( lind$sl ==x)  )

source("R/setLineageColors.r")
clColors <- seq(1, length(lineage))
clColors[ colorInd] <- setLineageColors( lineage, lineageId, MC_linId, IP_linId)

vln <- VlnPlot(seuratObj, gene, ident.include = slingObj@lineages[[lineageId]], do.return = TRUE)

png( file.path( plotDir, "Vln", paste0( gene, "_l_", lineageId, "_vln.png")))
	plot(vln 
		+ scale_x_discrete(limits = slingObj@lineages[[lineageId]])
		+ scale_fill_manual(values = clColors)
		)
dev.off()

}
