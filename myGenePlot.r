myGenePlot <- function(GeneCellType, Gene_1, Gene_2, Types = "All", Cols = "black", png = FALSE){

require("Seurat")

#GeneCellType is a Seurat object. Cell types are in the @ident and cell names are in @cell.names slots.

if(length(Types) == 1){if(Types == "Alltypes"){Types <- GeneCellType@ident}} 

names(Cols)<-Types

if(length(Types) != length(Cols)){return("Colors and Types do not match")} 
Gene_1e <- FetchData(GeneCellType, Gene_1)
Gene_2e <- FetchData(GeneCellType, Gene_2)

if(png){
	PlotFileName <- file.path(getwd(), "Plot", paste(Gene_2, "_vs_", Gene_1, ".png", sep = ""))
	png(filename = PlotFileName)
	plot(Gene_2e ~ Gene_1e, 
		xlim = c(min(Gene_1e), max(Gene_1e)), 
		ylim = c(min(Gene_2e), max(Gene_2e)), 
		xlab = Gene_1, ylab = Gene_2, 
		main = paste(Gene_2, " vs ", Gene_1, sep = ""))


		for (type in Types){points(Gene_2e[GeneCellType@ident==type]~Gene_1e[GeneCellType@ident==type], col = Cols[type], pch = 19,  
					xlim = c(min(Gene_1e), max(Gene_1e)), ylim = c(min(Gene_2e), max(Gene_2e)))} #over Types
	dev.off()
	} # if png is TRUE
	else {
		plot(Gene_2e ~ Gene_1e, 
		xlim = c(min(Gene_1e), max(Gene_1e)), 
		ylim = c(min(Gene_2e), max(Gene_2e)), 
		xlab = Gene_1, ylab = Gene_2, 
		main = paste(Gene_2, " vs ", Gene_1, sep = ""))	

		for (type in Types){points(Gene_2e[GeneCellType@ident==type]~Gene_1e[GeneCellType@ident==type], col = Cols[type], pch = 19,  
			xlim = c(min(Gene_1e), max(Gene_1e)), ylim = c(min(Gene_2e), max(Gene_2e)))} #over Types
	} # if png is FALSE

}	#myGenePlot
