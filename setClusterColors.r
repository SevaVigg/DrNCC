setClusterColors <- function(seuratObj, clustType){

#This snippet sets a color code for all cells according to their clusters from obj@ident. Used to draw slingshot lineages

cellColors <- rep("grey", length(seuratObj@ident))

cellColors[ seuratObj@ident == clustType["I"]] <- "cyan"
cellColors[ seuratObj@ident == clustType["M"]] <- "black"
cellColors[ seuratObj@ident == clustType["Tl"]]<- "red"

setClusterColors <- cellColors
}
