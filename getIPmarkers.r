# This snippet uses file SOM codes 45 14 03 prepared by SOM clasterer to identify iridophore markers

require("data.table")

SomCodes <- read.table("Source/SOME_codes_45_14_03.txt", header = T)

Clusters <- read.table("Source/Clusters_gSOM_45_14_03.txt", header = T)

Clusters <- as.data.table(Clusters)
Clusters[, type := tstrsplit(name, "\\.")[1]]

MaxIP_cl <- Clusters[,.(N_IP = sum(type == "IP")), class][order(-N_IP)][1, class]

IPgenes <- unlist(sort(SomCodes[MaxIP_cl,], decreasing = T))
IPmarkers <- IPgenes[which(IPgenes > log10(2))]

OutFileName 	<- file.path(getwd(), "Source", "Markers", "SOM_IPmarkers.txt")
OutFile		<- file(OutFileName, open = "w")

write.table(IPmarkers, file = OutFile, col.names = F)

close(OutFile) 		


