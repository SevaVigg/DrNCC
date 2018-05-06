setInitColors <- function(seuratObj){

#This functions sets a color code for initial cells 
 
cellColors <- rep("grey", length(seuratObj@ident))

cellColors[ which(levels(seuratObj@ident) == "I")] = "cyan"
cellColors[ which(levels(seuratObj@ident) == "M")] = "black"
cellColors[ which(levels(seuratObj@ident) == "m6")] = "blue"
cellColors[ which(levels(seuratObj@ident) == "W2")] = "magenta"
cellColors[ which(levels(seuratObj@ident) == "Tl")] = "red"

setInitColors <- cellColors
}
