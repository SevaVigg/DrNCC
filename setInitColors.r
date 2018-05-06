setInitColors <- function(seuratObj){

#This snippet sets a color code for all initial cells 
 
cellColors <- rep("grey", length(seuratObj@ident))

cellColors[ grep("I",     names( seuratObj@ident))] = "cyan"
cellColors[ grep("M",     names( seuratObj@ident))] = "black"
cellColors[ grep("m618",  names( seuratObj@ident))] = "blue"
cellColors[ grep("mitfa", names( seuratObj@ident))] = "magenta"
cellColors[ grep("tails", names( seuratObj@ident))] = "red"

setInitColors <- cellColors
}
