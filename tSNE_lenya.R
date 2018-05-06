#loading table into list
count <- 0
getTable <- function (filename) {
	data2 <- read.csv(filename, sep=",", header=TRUE, stringsAsFactors=FALSE, check.names = FALSE)
#	print(colnames(data2))
	if (!is.null(data2[["Gene.Name"]])) {
		data2 <- data2[nchar(colnames(data2)) != 0]
#		print(strsplit(colnames(data2), "_"))
		data2 <- data2[!data2[,"Gene.Name"] == "",]
		rownames(data2)<-data2[,"Gene.Name"]
		cols <- paste(colnames(data2), as.character(count), sep="_")
		cols <- gsub(" new ", " ", cols)
		cols <- gsub("hpf", "HPF", cols)
		cols <- gsub(" mixed _", "", cols)
		cols <- gsub("HPF_", "HPF", cols)
		cols <- gsub("30HPF13", "30HPF 13", cols)
		cols <- gsub("21 HPF", "21HPF", cols)
		cols <- gsub("study24HPF", "study 24HPF", cols) 
		cols <- gsub("study48HPF07", "study 48HPF 07", cols) 
		cols <- gsub("study30HPF", "study 30HPF", cols) 
		colnames(data2)<-cols
	}
	else
	{
		data2 <- data2[nchar(colnames(data2)) != 0]
		cols<-paste(colnames(data2), as.character(count), sep="_")
		rows <- data2[nchar(data2[,1]) != 0,1]
		data2<-data2[nchar(data2[,1]) != 0,]
#		print(length(rownames(data2)))
#		print(rows)
		rownames(data2)<-rows
#		print(cols)
#		cols <- 
		cols <- gsub(" new ", " ", cols)
		cols <- gsub("hpf", "HPF", cols)
		cols <- gsub(" mixed _", "", cols)
		cols <- gsub("HPF_", "HPF ", cols)
		cols <- gsub("30HPF13", "30HPF 13", cols)
		cols <- gsub("21 HPF", "21HPF", cols)
		cols <- gsub("study24HPF", "study 24HPF", cols) 
		cols <- gsub("study48HPF", "study 48HPF", cols) 
		cols <- gsub("study30HPF", "study 30HPF", cols) 
		print("!!!!!!!")
		print(data2[1])
		print("!!!!!!!")
#		cols <- lapply(cols, function (x) { return(toupper(x)) })
#		print(cols)
		colnames(data2)<-cols
	#	colnames(data2)<-paste(colnames(data2), as.character(count), sep="_")
	}
	count <<- count+1
	print(count)
	return(data2[,grepl("study",names(data2[1,]))])
}
filenames <- list.files("Zerbra_fish_2", pattern="study(.*)csv$", full.names=TRUE)
ldf <- list()
ldf <- c(ldf,lapply(filenames, getTable))
z<-as.data.frame(ldf)
library(plyr)
m <- laply(z, function (x) { return (as.double(x)) })
colnames(m)<-rownames(ldf[[1]])
n<-strsplit(paste(laply(ldf, function(x) { return(paste(colnames(x), collapse='!!')) } ), collapse = '!!'), "!!")
rownames(m) <- n[[1]]
library(Seurat)
groups<-laply(strsplit(rownames(m), " "), function(x) { return(x[3])} )
print(groups)
library(limma)
y2 <- removeBatchEffect(t(log(m)), groups)
zebrafish <- new("seurat", raw.data = y2)
zebrafish <- Setup(zebrafish, project="Zebrafish", min.genes=60, names.delim=" ",  names.field=3)
zebrafish <- MeanVarPlot(zebrafish, fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = -12, x.high.cutoff = 12, y.cutoff = 0.1, do.contour = F, do.plot=T)
zebrafish <- PCA(zebrafish)
zebrafish <- JackStraw(zebrafish, num.replicate = 100, do.print = FALSE)
zebrafish <- RunTSNE(zebrafish, dims.use = 1:5, add.iter=1000)
