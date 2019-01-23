
if (!require("vsn")){
if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
	BiocManager::install("vsn", version = "3.8")
}

if(!require("NanoStringNorm")){
  install.packages("NanoStringNorm", dependencies = TRUE)
}

require(preprocessCore)
require(data.table)

plotDir <- file.path(getwd(), "Plot")
resDir	<- file.path(getwd(), "Res")

qualContLogFile	<- file(paste0(resDir, .Platform$file.sep, "qualContLog.txt"), open = "w")

Genes 	<- read.csv( file = paste0("Res", .Platform$file.sep, "expressionTableDedup.csv") , stringsAsFactors = FALSE)
Cells	<- read.csv( file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedup.csv") , stringsAsFactors = FALSE)
Probes	<- read.csv( file = paste0("Res", .Platform$file.sep, "ProbesDescripitonsDedup.csv") , stringsAsFactors = FALSE)

geneNames <- Genes[,1]
cellNames <- Cells[,1]

rownames(Cells)	<- cellNames
rownames(Genes)	<- geneNames
Cells[,1]	<- NULL
Genes[,1]	<- NULL

# We have a number of x18 cells with missing mitfa values; these values are updated as medians
# X18_mitfa 	<- Genes["mitfa", grep("X18", colnames(Genes))]
# X18_mitfa_med	<- median(as.numeric(X18_mitfa["mitfa", !is.na(X18_mitfa[1,])]))
# Genes["mitfa", is.na(Genes["mitfa",])] <- X18_mitfa_med
# cat(file = qualContLogFile, ncol(X18_mitfa), " cells with missing mitfa values updated as medians\n")

#remove mitfa gene

colnames(Genes) <- colnames(Cells) <- paste0(Cells["hpf",],"_", Cells["CellType",])

geneMatrix	<- as.matrix(Genes)

normGeneMatrix	<- normalize.quantiles(geneMatrix)

rownames(normGeneMatrix) <- Probes[, "Gene.Name"]
colnames(normGeneMatrix) <- colnames(Cells)


rpl13_Distr	<- log( unlist( normGeneMatrix["rpl13", ]))

png( file.path( plotDir, "rpl13_Distr_plot.png"))
	plot(sort(rpl13_Distr))
	abline( h = quantile( rpl13_Distr, 0.05), col = "red")
	abline( h = quantile( rpl13_Distr, 1), col = "red")
	abline( h = min( rpl13_Distr[ grep("MC", names(rpl13_Distr))]), col = "black")
	abline( h = max( rpl13_Distr[ grep("MC", names(rpl13_Distr))]), col = "black")
	abline( h = max( rpl13_Distr[ grep("IP", names(rpl13_Distr))]), col = "cyan")
	abline( h = min( rpl13_Distr[ grep("IP", names(rpl13_Distr))]), col = "cyan")
dev.off()


rpl13_keep	<- which( rpl13_Distr > quantile(rpl13_Distr, 0.01) )

Genes_f			<- Genes[,rpl13_keep]
Cells_f			<- Cells[,rpl13_keep]
geneMatrix_f		<- as.matrix(Genes_f)
normGeneMatrix_f	<- normalize.quantiles( geneMatrix_f)

rownames(normGeneMatrix_f) <- Probes[, "Gene.Name"]
colnames(normGeneMatrix_f) <- colnames(Cells_f)


cat(file = qualContLogFile, ncol(Genes_f), " cells remaining after removing cells with rpl13 dropouts\n")

Kanamycin_Distr	<- log( unlist( normGeneMatrix_f["Kanamycin Pos", ]))

png( file.path( plotDir, "Kanamycin_Distr_plot.png"))
	plot(sort( Kanamycin_Distr))
	abline( h = quantile( Kanamycin_Distr, 0.01), col = "red")
	abline( h = quantile( Kanamycin_Distr, 0.95), col = "red")
	abline( h = min( Kanamycin_Distr[ grep("MC", names( Kanamycin_Distr))]), col = "black")
	abline( h = max( Kanamycin_Distr[ grep("MC", names( Kanamycin_Distr))]), col = "black")
	abline( h = max( Kanamycin_Distr[ grep("IP", names( Kanamycin_Distr))]), col = "cyan")
	abline( h = min( Kanamycin_Distr[ grep("IP", names( Kanamycin_Distr))]), col = "cyan")
dev.off()

KanamycinPos_keep	<- which( Kanamycin_Distr > quantile(Kanamycin_Distr, 0.01) )

Genes_ff		<- Genes_f[,KanamycinPos_keep]
Cells_ff		<- Cells_f[,KanamycinPos_keep]
geneMatrix_ff		<- as.matrix(Genes_ff)
normGeneMatrix_ff	<- normalize.quantiles( geneMatrix_ff)

rownames(normGeneMatrix_ff) <- Probes[, "Gene.Name"]
colnames(normGeneMatrix_ff) <- colnames(Cells_ff)


cat(file = qualContLogFile, ncol(Genes_ff), " cells remaining after removing cells with Kanamycin Pos dropouts\n")

batches		<- unique(unlist(Cells_ff["FileName",]))
batchDates	<- unlist(lapply(batches, function(x) Cells_ff["dateEx", grep(x, Cells_ff["FileName",], fixed = TRUE)][1]))
batchTypes	<- unlist(lapply(batches, function(x) Cells_ff["CellType", grep(x, Cells_ff["FileName",], fixed = TRUE)][1]))
batchHpf	<- unlist(lapply(batches, function(x) Cells_ff["hpf", grep(x, Cells_ff["FileName",], fixed = TRUE)][1]))

batchVals	<- lapply(batches, function(x) geneMatrix[, which(Cells_ff["FileName", ]==x)])
batchNormVals 	<- lapply(batches, function(x) normGeneMatrix[, which(Cells_ff["FileName", ]==x)])

names(batchNormVals) <- as.character(batches)
names(batchVals)     <- as.character(batches)

batchMedVals	<- unlist(lapply(batchVals, median))
batchMedNormVals<- unlist(lapply(batchNormVals, median))

qualMatrix	<- as.data.frame(rbind(batchHpf, batchDates, batchTypes, batchMedVals, batchMedNormVals), stringsAsFactors = FALSE)
names(qualMatrix)	<- batches
qualMatrix	<- qualMatrix[, order(as.numeric(qualMatrix["batchMedNormVals",]))]

png(paste0(plotDir, .Platform$file.sep,"batchMedians.png"))
	plot(as.numeric(qualMatrix["batchMedNormVals", ]), main = "Medians of batches", xlab = "Batches", ylab = "CountMedian")
	abline(h = 40, col = "red")
	abline(h = min(as.numeric(qualMatrix["batchMedNormVals", qualMatrix["batchTypes",]=="MC"])), col = "blue")
dev.off()

batchVecs <- lapply(batchVals, function(x) log(as.vector(x)))

png(paste0(plotDir, .Platform$file.sep, "LogNotNormedBatchBoxPlot.png"), width = 960)
	boxplot(batchVecs[colnames(qualMatrix)], 
		main = "Not normalized, not filtered batches (log-scale)", 
		las =2 , 
		cex.axis = 0.7, 
		boxwex = 1.2, 
		at = seq(1, 2*length(batches), 2))
dev.off()

batchProbl 	<- qualMatrix[, which(as.numeric(qualMatrix["batchMedNormVals",])<40)]	#threshold to keep MC (MedNormVal ~ 60)
cellsProbl_Ind	<- which(Cells_ff["FileName",] %in% colnames(batchProbl))


Genes_fff		<- Genes_ff[,-cellsProbl_Ind]
Cells_fff		<- Cells_ff[,-cellsProbl_Ind]

cat(file = qualContLogFile, ncol(Genes_fff), " cells remaining after removing batches with low medians\n")

posSpikes	<- c(128, 32, 8, 2, 0.5, 0.125)						# Pos probes in fM
posProbes	<- c("POS_A", "POS_B", "POS_C", "POS_D", "POS_E", "POS_F")

coefs		<- apply( log2(Genes_fff[posProbes,]), 2, function(x) lm( x~log2(posSpikes))$coefficients[2] )

sortLogCoefs	<- sort(coefs)

png(paste0(plotDir, .Platform$file.sep, "LogSortPosCoefs.png"))
	plot(sortLogCoefs, main = "Positive quality control coef values")
	abline( h =1, col = "green")
dev.off()

keepCoefs <- which(coefs > 0.9 & coefs < 1.1)

Genes_ffff	<- Genes_fff[,keepCoefs]
Cells_ffff	<- Cells_fff[,keepCoefs]

cat(file = qualContLogFile, ncol(Genes_ffff), " cells remaining after removing cells with poor positive control values\n")

normIteration <- 0
repeat{									       #iterations over background level								      
normIteration <- normIteration + 1
cat("Starting iteration ", normIteration, "\n")

NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_ffff, stringsAsFactors = FALSE)
colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")

NanoTable[c("Kanamycin Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
NanoTableNormed <- NanoStringNorm(x = NanoTable, CodeCount = "sum", Background = "mean", SampleContent = "housekeeping.geo.mean")

nf	 	<- NanoTableNormed$sample.summary.stats.norm$pos.norm.factor	#recommended values for normfactors are between 0.3 and 3
bg	 	<- NanoTableNormed$sample.summary.stats.norm$background.level	#recommended values for background are within 3 sd from the mean
sm	 	<- NanoTableNormed$sample.summary.stats.norm$Sample.Missing	#recommended values for the number of missing samples is less than 0.9
rc	 	<- NanoTableNormed$sample.summary.stats.norm$rna.content	#recommended values for RNA content are within 3 sd from the mean
nf_drop		<- which( nf < 0.3 | nf > 3)
sm_drop		<- which( sm > 0.9 )
bg_drop		<- which( abs( bg - mean(bg)) > 3*sd(bg))
rc_drop		<- which( abs( rc - mean(rc)) > 3*sd(rc))

cat( file = qualContLogFile, "Iteration ", normIteration, length(nf_drop), " cells with norm factors out of the recommended band\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(bg_drop), " cells with the backgrounds out of the recommend band\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(sm_drop), " cells with the number of missing samples less than 0.9\n")
cat( file = qualContLogFile, "Iteration ", normIteration, length(rc_drop), " cells with RNA content deviating too much from the mean\n")

"%u%"		<- union							#union() supports only two arguments, define intersect operation
norm_drop	<- nf_drop %u% sm_drop %u% bg_drop %u% rc_drop

if( length(norm_drop) == 0){break}

Genes_ffff	<- Genes_ffff[, -norm_drop]
Cells_ffff	<- Cells_ffff[, -norm_drop]

}										# repeat

cat( file = qualContLogFile, ncol(Genes_ffff), " cells remaining after removing cells with poor normalization statistics\n")

NanoTable	<- cbind(Probes[,"Class.Name"], Probes[,"Gene.Name"], Probes[,"Accession.."], Genes_ffff, stringsAsFactors = FALSE)
colnames(NanoTable)[1:3] <- c("Code.Class", "Name", "Accession")
NanoTable[c("Kanamycin Pos", "rpl13"), "Code.Class"] <- "Housekeeping"
Genes_n		 <- NanoStringNorm(x = NanoTable, CodeCount = "sum", Background = "mean", SampleContent = "housekeeping.geo.mean", return.matrix.of.endogenous.probes = TRUE)

geneMatrix_n	<- as.matrix(Genes_n)

batches_n	<- unique(unlist(Cells_ffff["FileName",]))
batchDates_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["dateEx", grep(x, Cells_ffff["FileName",], fixed = TRUE)][1]))
batchTypes_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["CellType", grep(x, Cells_ffff["FileName",], fixed = TRUE)][1]))
batchHpf_n	<- unlist(lapply(batches_n, function(x) Cells_ffff["hpf", grep(x, Cells_ffff["FileName",], fixed = TRUE)][1]))

batchVals_n	<- lapply(batches_n, function(x) geneMatrix_n[, which(Cells_ffff["FileName", ]==x)])

names(batchVals_n) <- as.character(batches_n)

batchMedVals_n	<- unlist(lapply(batchVals_n, median))

qualMatrix_n	<- as.data.frame(rbind(batchHpf_n, batchDates_n, batchTypes_n, batchMedVals_n), stringsAsFactors = FALSE)
names(qualMatrix_n)	<- batches_n
qualMatrix_n	<- qualMatrix_n[, order(as.numeric(qualMatrix_n["batchMedVals_n",]))]

png(paste0(plotDir, .Platform$file.sep,"NormalizedBatchMedians.png"))
	plot(as.numeric(qualMatrix_n["batchMedVals_n", ]), main = "Medians of normalized batches", xlab = "Batches", ylab = "CountMedian")
dev.off()

batchVecs_n <- lapply(batchVals_n, function(x) log(as.vector(x)+1))

png(paste0(plotDir, .Platform$file.sep, "LogNormalizedBatchBoxPlot.png"), width = 960)
	boxplot(batchVecs_n[colnames(qualMatrix_n)], 
		main = "Normalized and filtered batches (log-scale)", 
		las =2 , 
		cex.axis = 0.7, 
		boxwex = 1.2, 
		at = seq(1, 2*length(batches_n), 2))
dev.off()


write.table(Genes_n, file = paste0(resDir, .Platform$file.sep, "NormalizedExTable.csv"), sep = "\t")



cells_dt 	<- as.data.table( t(Cells_ffff))
cells_tbl	<- cells_dt[, .N, .(CellType,hpf)]

write.table(cells_tbl, file=qualContLogFile, sep = "\t")

close(qualContLogFile)	

write.table( Genes_ffff, file = paste0("Res", .Platform$file.sep, "expressionTableDedupQC.csv"), sep = "\t" )
write.table( Cells_ffff, file = paste0("Res", .Platform$file.sep, "cellDescripitonsDedupQC.csv"),sep = "\t" )









 


		
 

