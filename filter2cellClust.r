filter2cellClust <- function(object, clustArray){

#this function removes cluster from clustArray (containing ident nums) that contain only two cells. This is the patch for VlnPlot, which incorrectly
#displays colors of clusters that contain only two cells

require(Seurat)

cellNums <- sapply( clustArray, function(x) length( WhichCells( object, x)))
resArray <- clustArray[ which( cellNums >2) ]
return( resArray)
}
