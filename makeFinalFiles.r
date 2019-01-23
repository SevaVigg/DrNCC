#Preparing clean  data for WT, WTandSox10, WTandAllmutants

cleanDataDir	<- file.path(getwd(), "cleanData")
dir.create( cleanDataDir, showWarnings = FALSE)



if(!require(DrImpute)){
  library(devtools)
  install_github('gongx030/DrImpute')		
  library("DrImpute")
}

imputeDir	<- file.path(getwd(), "Imputed")
dir.create( imputeDir, showWarnings = FALSE)

imputeWTDir 	<- file.path( imputeDir, "WT")
dir.create( imputeWTDir, showWarnings = FALSE)

imputeWTResDir <- file.path( imputeWTDir, "Res")
dir.create( imputeWTResDir, showWarnings = FALSE)

dirQCres <- file.path( getwd(), "QualityControl", "Res") 

genesQC	<- read.table( file = file.path( dirQCres, "NormalizedExTable.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )
cellsQC	<- read.table( file = file.path( dirQCres, "cellDescripitonsDedupQC.csv"), sep = "\t", stringsAsFactors = FALSE, check.names=FALSE )

rownames(genesQC)[which(rownames(genesQC)=="Kanamycin Pos")] <- "Kanamycin_Pos"

allGenes	<- rownames(genesQC)
logExps 	<- log2(1+genesQC)

#drop mutants
logExpsWT	<- logExps[, setdiff( seq_along(colnames(logExps)), grep("(sox10|mitfa)", colnames(logExps)))]
cellsImpWT	<- cellsQC[, setdiff( seq_along(colnames(Cells)), grep("(sox10|mitfa)", colnames(Cells)))]

#and now impute for dropouts

seedDir	<- file.path( imputeWTResDir, "seed")
dir.create( seedDir, showWarnings = FALSE)

randSeed <- as.integer(Sys.time())
set.seed( randSeed )
cat( file = file.path(seedDir, "seed.txt"), randSeed)

logExpsImpWT <- DrImpute(logExpsWT)

write.table( logExpsImpWT, file = file.path( imputeWTResDir, "logExpTableDedupImpWT.csv"), sep = "\t" )
write.table( cellsImpWT, file = file.path( imputeWTResDir, "cellDescripitonsImpWT.csv"), sep = "\t" )
