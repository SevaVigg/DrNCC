

#This snippet requires MIGN seurat object

plotDir <- file.path(getwd(), "Plot")
plotVlnDir <- file.path(getwd(), "Plot", "Vln")
plotVlnunNormDir <- file.path(getwd(), "Plot", "Vln", "unNormalized")

resDir	<- file.path(getwd(), "Res")

dir.create( plotDir, showWarnings = FALSE)
dir.create( plotVlnDir, showWarnings = FALSE)
dir.create( plotVlnunNormDir, showWarnings = FALSE)

drawDir <- plotVlnunNormDir

lapply(rownames( MIGN@data), function(x){
	x <- unlist(x)
	MainText <- paste0(x, " unNormed")
	png(filename = paste0(drawDir, .Platform$file.sep, x, "_vln.png"))
	p <- VlnPlot(MIGN, x, cols.use = cellColors)
	plot(p)
	dev.off()
	})
