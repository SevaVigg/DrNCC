#This snippet assignes cell colors accordingly to their cell types
setCellTypeColors <- function(seuratObj){


if(!require("colorspace")){
install.packages("colorspace")
library(colorspace)}


#basic colors are shades of green, according to the clustering ident

numColors		<- length( levels( seuratObj@ident))
cellColors		<- integer( numColors)
names(cellColors)	<- levels( seuratObj@ident) 
standColors	<- grep( "[0-9][0-9]", names(cellColors), value = TRUE)
cellColors[standColors] <- sequential_hcl( length( standColors) , "BluGrn", alpha = 0.5, rev = TRUE)

names(cellColors) <- levels( unique(seuratObj@ident))

cellColors[ "I" ]   = "cyan"
cellColors[ "M" ]   = "black"
cellColors[ "m6" ]  = "gold"
cellColors[ "W2" ]  = "dodgerblue"
cellColors[ "Tl" ]  = "red"

#make tails red according to their "snail2" expression

return(cellColors)

}
