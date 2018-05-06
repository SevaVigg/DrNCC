plotDir	<- file.path(getwd(), "Plot", "MIGNVlnPlots")
dir.create( plotDir, showWarnings = FALSE)


genes <- c("mitfa", "ltk", "phox2b", "snail2", "sox9b", "pax7a", "pax7b", "mbpa", "sox10", "mlphb", "pnp4a", "tfec", "tyrp1b")

for (igene in genes){
	cat(igene, "\n")
	png( file.path(plotDir, paste0(igene, ".png")))
 		a <- VlnPlot(MIGN, cols.use = colorClusters(MIGN)[as.numeric(filter2cellClust(MIGN, levels(MIGN@ident)))], features.plot = igene)
		plot(a)	
 	dev.off()
 }


