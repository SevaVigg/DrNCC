colorClusters 		<- function(seurObj){
	clCols 		<- rep("grey", length(levels(seurObj@ident)))
	clCols[1]  	<- "black"
	clCols[2]	<- "blue"
	clCols[4:6]	<- "red"
	clCols[7:10]	<- "green"
	clCols[11:17]	<- "orange"
	clCols[18:25]	<- "grey"
return(clCols)
}
