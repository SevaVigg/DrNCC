allGenesVlns <- function( seuratObj, ident, genes,  titleStr){

if( !require("data.table"))install.packages("data.table")

cells 	<- WhichCells( seuratObj, ident)
texps	<- as.data.table( t( logExps[ genes, cells]))
dtExps	<- melt(texps)
names(dtExps) <- c("gene", "expression")

ggplot(dtExps, aes( x = gene, y = expression, fill = factor(gene))) + geom_violin(scale = "width", adjust = 0.5, show.legend = FALSE) + theme(axis.text.x = element_text(angle = 90, size = 20), axis.text.y = element_text(size = 20)) +labs( title = titleStr)->p
plot(p)
}
