setClusterColor <- function(seuratObj, clustType, Lineage = NULL){

#This snippet sets a color code for all cells according to their clusters from obj@ident. Used to draw slingshot lineages

clusters	<- unique(seuratObj@ident)

cellColors 	<- rep("grey", length( seuratObj@ident)) # initiate all cells as grey

#greyClust	<- grep("[0-9].", clusters)		 # standard cluster types
#colClust	<- setdiff( clusters, greyClust)	 # special cluster types

#greyPalette	<- grey.colors( length(greyClust))

#for (x in greyClust) cellColors[ clusters == x] <- greyPalette[ which(greyClust == x)]

cellColors[ seuratObj@ident == clustType["I"]] <- "cyan"
cellColors[ seuratObj@ident == clustType["M"]] <- "black"
cellColors[ seuratObj@ident == clustType["Tl"]]<- "red"
cellColors[ seuratObj@ident == clustType["m6"]]<- "pink"
cellColors[ seuratObj@ident == clustType["W2"]]<- "magenta"
if ( !is.null(Lineage)) 
	{cellColors[ seuratObj@ident == colnames(Lineage)[ncol(Lineage) - 1]] <- "#FF99FF"}

setClusterColor <- cellColors
}
