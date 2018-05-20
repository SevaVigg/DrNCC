library(leaflet)
library(ggplot2)
for (gene_name in rownames(logExps)) {
print(paste(gene_name, ".png", sep=""))
spline_x <- as.double(curves$curve29)
spline_y <- as.double(logExps[gene_name,names(curves$curve29)])
colors <- unlist(lapply(X=strsplit(x=names(curves$curve29), split="_"), FUN=function(x) { 
t1 <- strsplit(x=x[2], split="[.]")
if (t1[[1]][1] == "MC")
	return("MC")
if (t1[[1]][1] == "tails")
	return("tails")
if (t1[[1]][1] == "IP")
	return("IP")
if (t1[[1]][1] == "m618")
	return("m618")
if (t1[[1]][1] == "mitfa w2")
	return("mitfa w2")
return(x[1])}))
CellType <- colorNumeric(c("blue", "green"), 1:12)(as.integer(factor(colors)))#factor(colors)
CellType[colors == "MC"] <- "black"
CellType[colors == "IP"] <- "cyan"
CellType[colors == "m618"] <- "magenta"
CellType[colors == "tails"] <- "red"
CellType[colors == "mitfa w2"] <- "pink"
names(CellType) <- colors
p <- ggplot(data.frame(x=spline_x, y=spline_y, color=colors), aes(x, y))
p <- p+geom_point(aes(colour = colors), size = 3)+ ggtitle(paste(gene_name, " gene"))+geom_smooth(method="loess", span = 0.3)+theme(legend.title=element_blank(), plot.title = element_text(size = 40, face = "bold"), axis.title=element_text(size=20,face="bold"))+xlab("Pseudotime")+ylab("LogExpression")+scale_color_manual(values=CellType)
png(paste(gene_name, ".png", sep=""), height=1000, width=1000)
print(p)
dev.off()
}
