#####################################################################################################################
# This function gets a seurat object 'object' and runs SNN.r function of ranked SNN tree created by Chen Xu
# and downloaded from http://bioinfo.uncc.edu/SNNCliq/#download
#
# Written by Vsevolod Makeev at VIGG, started 10.04.2017, finished ...
#
#####################################################################################################################
require("Seurat")
require("Matrix")

source(file.path(getwd(),"R","SNN.r"))

buildSNNcliq <- function( object, k, distance = "euclidian"){

if (missing(k)){k = 3}

data_matrix <- t( as.matrix( object@data))
ncells <- length( object@cell.names)

tc <- textConnection( "SNN_out", open = "w")		     #to intersept output of SNN.r, which is "write.table"	
	SNN( data_matrix, tc, k, distance = distance) #Runs SNN.r , now tc contains a list, "\t" separated
close(tc)

SNN_list <- transpose( strsplit( SNN_out, "\t"))		      #list of three columns. As character
SNN_sparse <- sparseMatrix( i = as.integer( SNN_list[[1]]), j = as.integer( SNN_list[[2]]), x = as.numeric( SNN_list[[3]]), dims = c(ncells,ncells))
SNN_dense  <- as.matrix( SNN_sparse)
SNN_dense  <- SNN_dense + t( SNN_dense)

row.names( SNN_dense) <- colnames( SNN_dense) <- object@cell.names

object@snn.dense <- SNN_dense

return(object)
} 
