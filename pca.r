Ded<-read.table("Source/Deduplicated.csv", header =T, sep =",")
row.names( Ded)<-Ded[,1]
Ded<-Ded[,-1]


types<-c("X18", "X21",  "X24",  "X30",  "X36",  "X48",  "X60",  "X72",  "X30W", "I", "M","silva_neg_mlphb_pos")



plot( hclust( dist( log( Ded+1))))

tDed<-as.data.frame( t( Ded))
logpca<-prcomp( log( tDed+1), center = T, scale =T)

pc2<-logpca$x[,"PC2"]
pc3<-logpca$x[,"PC3"]
pc4<-logpca$x[,"PC4"]

for( type in types) {
type_pc2<-logpca$x[grep( type, rownames( logpca$x)), "PC2"]
type_pc3<-logpca$x[grep( type, rownames( logpca$x)), "PC3"]
type_pc4<-logpca$x[grep( type, rownames( logpca$x)), "PC4"]

m_pc2<-logpca$x[grep( "M", rownames( logpca$x)), "PC2"]
i_pc2<-logpca$x[grep( "I", rownames( logpca$x)), "PC2"]
m_pc3<-logpca$x[grep( "M", rownames( logpca$x)), "PC3"]
i_pc3<-logpca$x[grep( "I", rownames( logpca$x)), "PC3"]
m_pc4<-logpca$x[grep( "M", rownames( logpca$x)), "PC4"]
i_pc4<-logpca$x[grep( "I", rownames( logpca$x)), "PC4"]

silva_neg_mlphb_pos_pc2<-logpca$x[tDed[tDed$"silva"<200 & tDed$"mlphb">5000,],"PC2"]
silva_neg_mlphb_pos_pc3<-logpca$x[tDed[tDed$"silva"<200 & tDed$"mlphb">5000,],"PC3"]
silva_neg_mlphb_pos_pc4<-logpca$x[tDed[tDed$"silva"<200 & tDed$"mlphb">5000,],"PC4"]


plotFileName<-file.path( getwd(), "Plot", paste(type, "_PC2_PC3.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC2"], logpca$x[,"PC3"], xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), main = paste(type, "_PC2_PC3", sep = ""))
	points(type_pc2, type_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "green", pch = 19)
	points(i_pc2, i_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "red", pch = 19)
	points(m_pc2, m_pc3, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc3), max(pc3)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( getwd(), "Plot", paste(type, "_PC2_PC4.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC2"], logpca$x[,"PC4"], xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), main = paste(type, "_PC2_PC4", sep = ""))
	points(type_pc2, type_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "green", pch = 19)
	points(i_pc2, i_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "red", pch = 19)
	points(m_pc2, m_pc4, xlim = c(min(pc2), max(pc2)), ylim = c(min(pc4), max(pc4)), col = "blue", pch = 19)
dev.off()

plotFileName<-file.path( getwd(), "Plot", paste(type, "_PC3_PC4.png" ,sep=""))
png(filename = plotFileName)
	plot(logpca$x[,"PC3"], logpca$x[,"PC4"], xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), main = paste(type, "_PC3_PC4", sep = ""))
	points(type_pc3, type_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "green", pch = 19)
	points(i_pc3, i_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "red", pch = 19)
	points(m_pc3, m_pc4, xlim = c(min(pc3), max(pc3)), ylim = c(min(pc4), max(pc4)), col = "blue", pch = 19)
dev.off()

	} #plot over cell types

png("Plot/bipot_PC2_PC4.png")
	biplot(logpca, choices = c(2,4), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png("Plot/bipot_PC2_PC3.png")
	biplot(logpca, choices = c(2,3), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png("Plot/bipot_PC3_PC4.png")
	biplot(logpca, choices = c(3,4), cex = c(1,.5), xlabs = rep(".", nrow(logpca$x)) )
dev.off()

png("Plot/dendro.png")
	plot(hclust(dist(log(Ded+1))), cex = 0.7)
dev.off()

