plotDir <- file.path(getwd(), "Plot")
plotPCADir <- file.path(getwd(), "Plot", "PCA")
plotPCAunNormDir <- file.path(getwd(), "Plot", "PCA", "unNormalized")

resDir	<- file.path(getwd(), "Res")

dir.create( plotDir, showWarnings = FALSE)
dir.create( plotPCADir, showWarnings = FALSE)
dir.create( plotPCAunNormDir, showWarnings = FALSE)

drawDir <- plotPCAunNormDir

Genes 	<- read.table( file = paste0("Res", .Platform$file.sep, "expressionTableDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
Cells	<- read.table( file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)
Probes	<- read.table( file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE)

cells_ind 	<- order(as.numeric(Cells["hpf",]))
Genes		<- Genes[, cells_ind]
Cells		<- Cells[, cells_ind]

types 	<- unique(paste0(Cells["hpf",], "_", Cells["CellType",]))
types30 <-  c("X18", "X21",  "X24", "tails", "mitfa", "m618", "IP", "MC")

newtypes   <-  c("18", "21", "24", "Tl", "30", "W2", "m6", "36", "48", "I", "M", "60", "72")
newtypes30 <-  c("18", "21", "24", "Tl", "W2", "So", "m6", "I", "M")

cellColors <- c("lightgreen", "green", "darkgreen", "turquoise", "forestgreen", "orange", "red", "grey", "grey", "magenta", "blue", "grey", "grey")

png( paste0( drawDir, .Platform$file.sep, "GeneProbeDendroNotNormalized.png"))
	par(cex = 0.7)
	plot( hclust( dist( log (Genes +1))), main = "Cluster Dendrograme: Genes and NanoString Probes; No normalization", xlab = "Complete linkage")
dev.off()

tGenes<-as.data.frame( t( Genes))
logpca<-prcomp( log( tGenes+1), center = T, scale =T)

pc1	<- logpca$x[, "PC1"]
pc2	<- logpca$x[, "PC2"]
pc3	<- logpca$x[, "PC3"]
pc4	<- logpca$x[, "PC4"]

for( type in types) {

type_pc1	<-logpca$x[grep( type, rownames( logpca$x)), "PC1"]	
type_pc2	<-logpca$x[grep( type, rownames( logpca$x)), "PC2"]
type_pc3	<-logpca$x[grep( type, rownames( logpca$x)), "PC3"]
type_pc4	<-logpca$x[grep( type, rownames( logpca$x)), "PC4"]

m_pc1<-logpca$x[grep( "M", rownames( logpca$x)), "PC1"]
i_pc1<-logpca$x[grep( "I", rownames( logpca$x)), "PC1"]
m_pc2<-logpca$x[grep( "M", rownames( logpca$x)), "PC2"]
i_pc2<-logpca$x[grep( "I", rownames( logpca$x)), "PC2"]
m_pc3<-logpca$x[grep( "M", rownames( logpca$x)), "PC3"]
i_pc3<-logpca$x[grep( "I", rownames( logpca$x)), "PC3"]
m_pc4<-logpca$x[grep( "M", rownames( logpca$x)), "PC4"]
i_pc4<-logpca$x[grep( "I", rownames( logpca$x)), "PC4"]

plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC1_PC2.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC1"], logpca$x[,"PC2"], xlim = c(min(pc1), max(pc1)), ylim = c(min(pc2), max(pc2)), main = paste(type, "_PC1_PC2", sep = ""))
	points(type_pc1, type_pc2, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc2), max(pc2)), col = "green", pch = 19)
	points(i_pc1, i_pc2, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc2), max(pc2)), col = "red", pch = 19)
	points(m_pc1, m_pc2, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc2), max(pc2)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC1_PC3.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC1"], logpca$x[,"PC3"], xlim = c(min(pc1), max(pc1)), ylim = c(min(pc3), max(pc3)), main = paste(type, "_PC1_PC3", sep = ""))
	points(type_pc1, type_pc3, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc3), max(pc3)), col = "green", pch = 19)
	points(i_pc1, i_pc3, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc3), max(pc3)), col = "red", pch = 19)
	points(m_pc1, m_pc3, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc3), max(pc3)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC1_PC4.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC1"], logpca$x[,"PC4"], xlim = c(min(pc1), max(pc1)), ylim = c(min(pc4), max(pc4)), main = paste(type, "_PC1_PC4", sep = ""))
	points(type_pc1, type_pc4, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc4), max(pc4)), col = "green", pch = 19)
	points(i_pc1, i_pc4, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc4), max(pc4)), col = "red", pch = 19)
	points(m_pc1, m_pc4, xlim = c(min(pc1), max(pc1)), ylim = c(min(pc4), max(pc4)), col = "blue", pch = 19)
dev.off()


plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC2_PC3.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC2"], logpca$x[,"PC3"], xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), main = paste(type, "_PC2_PC3", sep = ""))
	points(type_pc2, type_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "green", pch = 19)
	points(i_pc2, i_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "red", pch = 19)
	points(m_pc2, m_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC2_PC4.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC2"], logpca$x[,"PC4"], xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), main = paste(type, "_PC2_PC4", sep = ""))
	points(type_pc2, type_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "green", pch = 19)
	points(i_pc2, i_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "red", pch = 19)
	points(m_pc2, m_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( drawDir,  .Platform$file.sep, paste(type, "_PC3_PC4.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC3"], logpca$x[,"PC4"], xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), main = paste(type, "_PC3_PC4", sep = ""))
	points(type_pc3, type_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "green", pch = 19)
	points(i_pc3, i_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "red", pch = 19)
	points(m_pc3, m_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "blue", pch = 19)
dev.off()

	} #plot over cell types

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC1_PC2.png"))
	biplot(logpca, choices = c(1,2), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC1_PC3.png"))
	biplot(logpca, choices = c(1,3), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC1_PC4.png"))
	biplot(logpca, choices = c(1,4), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC2_PC4.png"))
	biplot(logpca, choices = c(2,4), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC2_PC3.png"))
	biplot(logpca, choices = c(2,3), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png(paste0( drawDir,  .Platform$file.sep, "biplot_PC3_PC4.png"))
	biplot(logpca, choices = c(3,4), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()


