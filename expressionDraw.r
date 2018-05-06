expressionDraw <- function(geneName, seurat_object) {

    allGenes	<- rownames(seurat_object@data)
    seurat_object@data <- round(seurat_object@data, digits = 2)	
    data 	<- seurat_object@data[geneName,]
    min_d 	<- min(data)
    max_d 	<- max(data)
    exprLevels 	<- seq(min_d, max_d, 0.01)
    
	
    seurat_object@ident <- as.factor(as.character( rep(0, length(seurat_object@ident))))
    levels(seurat_object@ident)	<- as.factor( as.character (exprLevels))
    seurat_object@ident 	<- as.factor( as.character( data ))

#   nLevels 			<- length( levels( seurat_object@ident ) )
    nLevels			<- length( exprLevels)
    names(seurat_object@ident)	<- colnames(data)

 
#    png(paste0( "Plot/Cruel/Distribs/", geneName, ".png"))
 
    TSNEPlot( seurat_object, no.legend = TRUE, colors.use = colorRampPalette(c("lightcyan", "blue"))(nLevels))
    
#TSNEPlot( seurat_object, no.legend = TRUE, colors.use = cm.colors( nLevels))	
	
#    dev.off()
#    data_c <- lapply(data, function(x) {
#        print(grey.colors(100, start = 0.0001, end=0.9999)[as.integer(100*((min_d-x)/(min_d-max_d)))+1]) 
#        return (grey.colors(100, start = 0.0001, end=0.9999)[as.integer(100*((min_d-x)/(min_d-max_d)))+1]) 
#        })
#    png(filename)
#    plot(x=ipmc_4pc@dr$tsne@cell.embeddings[,1], y=ipmc_4pc@dr$tsne@cell.embeddings[,2], col = unlist(data_c), xlab = "", ylab = "", xlim=c(-20, 20), ylim=c(-30,30), type="p", pch=16)
#    dev.off()

}
