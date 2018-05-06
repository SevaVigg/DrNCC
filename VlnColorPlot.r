VlnColorPlot <- function(object, ...){
require(Seurat)

cellNums <- sapply(levels(object@ident), function(x) length( WhichCells( object, x)))
return(cellNums)

}
