# This script identifies trajectories through clusters obtained by 
# cell subsampling 

sourceFileName 	<- file.path( getwd(), "Subsampling", "15_samples.txt")
subS		<- read.table( sourceFileName, row.names = 1, stringsAsFactors = TRUE )
nClust		<- apply(subS, 2, max)
cat(nClust)
