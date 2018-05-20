setLineageColors <- function( lineage, lineageId, MC_linId, IP_linId){

clustColors 	<- colorRampPalette(c("red", "blue"))(length(lineage))
clustColors[1]	<- "red"
if ( lineageId == MC_linId ) clustColors[length(lineage)] <- colorRampPalette( c("black", "black"))(1)
if ( lineageId == IP_linId ) clustColors[length(lineage)] <- colorRampPalette( c("cyan", "cyan"))(1)

return( clustColors)
}
