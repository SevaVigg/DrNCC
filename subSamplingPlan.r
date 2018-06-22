# This script identifies trajectories through clusters obtained by 
# cell subsampling 


sourceFileName 	<- file.path( getwd(), "Subsampling", "15_samples.txt")
subS		<- read.table( sourceFileName, row.names = 1, stringsAsFactors = TRUE )
nClust		<- apply(subS, 2, max)
cat(nClust, "\n")

source("R/seuratNorm.r")
    
origIdent       <- ipmc@ident

ipmc            <- SubsetData(ipmc, ident.remove = c("W2", "m6"))

comps		<- 	15

dir.create( file.path( resDir, Subsampling, paste0( "c", comps)), showWarnings = FALSE)
subSdir 	<- 	file.path( resDir, Subsampling, paste0( "c", comps))

# this is clustering for all cells
ipmc            <- RunPCA(ipmc, pc.genes = rownames( ipmc@data), pcs.compute = comps, do.print = FALSE) 
ipmc            <- RunTSNE(ipmc, dims.use = 1:comps, theta = 0, perplexity = 15)

source("R/plotInitTypesPcaTsne.r")
plotInitTypesPcaTsne( ipmc, subSdir)

resolDec 	<- 60

source("R/getClusterTypes.r")
clustTypes      <- getClusterTypes(ipmc)
