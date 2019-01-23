setClusterVarColor <- function(seuratObj, clustType, Lineage = NULL){

# This snippet sets a color code for all cells according to the variation of expression in their
# clusters from obj@ident. Used to draw slingshot lineages

clusters	<- unique(seuratObj@ident)

cellColors 	<- rep("grey", length( seuratObj@ident)) # initiate all cells as grey
nCells		<- ncol(seuratObj@data)

totVar		<- sum( apply( seuratObj@data, 1, function(x) sd(x)/mean(x)))
varLevels	<- seq(0, totVar)
nLevels		<- length(varLevels)
colors.use 	<- colorRampPalette(c("blue", "red"))(nLevels)


for (cluster in seq_along(clusters)){
	clData 	<- seuratObj@data[ , WhichCells( seuratObj, cluster)]
	nCells	<- ncol(clData)
	clVar	<- sum( apply( clData, 1, function(x) sd(x)/mean(x)))
	clLevel <- floor(clVar)
	cellColors[ seuratObj@ident == cluster] <- colors.use[clLevel]
	}
setClusterColors <- cellColors
}

