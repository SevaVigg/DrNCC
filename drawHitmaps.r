heatclus1<- clust[order(pseudotime(sds)[,1], na.last = NA)]
Genes_nh_heat <- ipmc_48@data[,order(pseudotime(sds)[,1], na.last = NA)]
set1Cols <- brewer.pal(8, "Set1")
png("Plot\\heatMap_1_6_5_4.version_4.png",width=5000,height=1000)
heatmap.2(log(as.matrix(Genes_nh_heat)+1),trace="none",Colv=NULL,ColSideColors=as.character(set1Cols[heatclus1]),cexCol = 1.5,cexRow = 1.5,margins = c(11, 6))
legend("topright",      
    legend = unique(heatclus1),
    col = unique(as.character(set1Cols[(heatclus1)])), 
    lty= 1,             
    lwd = 5,           
    cex=.7
    )
dev.off()


heatclus1<- clust[order(pseudotime(sds)[,2], na.last = NA)]
Genes_nh_heat <- ipmc_48@data[,order(pseudotime(sds)[,2], na.last = NA)]
set1Cols <- brewer.pal(8, "Set1")
png("Plot\\heatMap_1_6_2_3.version_4.png",width=5000,height=1000)
heatmap.2(log(as.matrix(Genes_nh_heat)+1),trace="none",Colv=NULL,ColSideColors=as.character(set1Cols[heatclus1]),cexCol = 1.5,cexRow = 1.5,margins = c(11, 6))
legend("topright",      
    legend = unique(heatclus1),
    col = unique(as.character(set1Cols[(heatclus1)])), 
    lty= 1,             
    lwd = 5,           
    cex=.7
    )
dev.off()
