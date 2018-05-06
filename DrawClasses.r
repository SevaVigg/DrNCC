require(data.table)

Prefix = "Fuzzy_gSOM"


Cl_table_FileName <- file.path( getwd(), "Source", paste("fuzzy_and_gsom.log.txt", sep = ""))
Cl_table_file <- file( Cl_table_FileName, open = "r")

Clust_raw <- read.table(Cl_table_file, header = TRUE)
close(Cl_table_file)

Clust_dt <- as.data.table( Clust_raw)

types<-c("X18", "X21",  "X24",  "X30",  "X36",  "X48",  "X60",  "X72",  "X30W", "I", "M")

Cl_stats <- Clust_dt[, .(.N, X = mean(coordX), Y = mean(coordY)), by = class][order(class)]

lapply(types, function(x) Clust_dt[grep( x, as.character( name)), type := x] )
Clust_dt[,type := as.factor(type)]

#Numbers
colors <- rev(terrain.colors(100))
PlotFileName <- file.path( getwd(), "Plot", "Tanyas_genes", paste(Prefix, "_", "N.png", sep = ""))
dir.create(file.path( getwd(), "Plot", "Tanyas_genes"), showWarnings = FALSE)

png( filename = PlotFileName)
	plot( Cl_stats[,X], Cl_stats[,Y], xlim = c( min( Cl_stats[,X]),max( Cl_stats[,X])), ylim = c( min(Cl_stats[,Y]),max( Cl_stats[,Y])), main = "N", cex =3)
	points( Cl_stats[,X], Cl_stats[,Y], col = colors[1+99*round( Cl_stats[,N]/sum(Cl_stats[,N]),digits = 2)], pch = 19, cex = 3)
	text( Cl_stats[,X], Cl_stats[,Y], labels = Cl_stats[,class])
dev.off()

#Distributions

colors <- rev(terrain.colors(100))
types<-levels( Clust_dt[,type])
lapply( types, function(x) {
		x_perc<-Clust_dt[order(type), .(X = mean(coordX), Y = mean(coordY), type = levels(type), 
			percentage = as.vector(table(type))/.N), by=class][type==x][order(class)]
		PlotFileName <- file.path( getwd(), "Plot", "Tanyas_genes", paste(Prefix, "_", x, "_distr.png" ,sep = ""))
		png(filename = PlotFileName)
			plot(Cl_stats[,X], Cl_stats[,Y], cex = 3, main = x, xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min(Cl_stats[,Y]),max(Cl_stats[,Y])))
			points(x_perc[,X], x_perc[,Y], xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min( Cl_stats[,Y]),max(Cl_stats[,Y])), cex = 3 , 
				col = colors[1+99*round(x_perc[,percentage],digits = 2)], pch = 19)
			text(Cl_stats[,X], Cl_stats[,Y], xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min(Cl_stats[,Y]),max(Cl_stats[,Y])), labels = Cl_stats[,class])
		dev.off()
			})#lapply

cells_FileName <- file.path( getwd(), "Source", paste("Deduplicated.csv", sep = ""))
cells_File <- file( cells_FileName, open = "r")

cells <- read.table( cells_File, header = TRUE, sep = ",")
close( cells_File)

genes <- levels(cells[,1])

cells <- as.data.table( cells)




#genes = c("mycl1a", "silva", "tyrp1b", "alx4b", "ednrba", "fgfr3_v2", "foxp4", "her9", "id2a", "impdh1b", "pnp4a", "sox5", "kita", "dpf3", "mlphb", "oca2","id2a")

#genes = c("alx4b", "dpf3", "ednrba", "foxd3", "id2a", "impdh1b", "kita", "ltk", "mbpa", "mitfa", "mlphb", "oca2", "pax7b", "pnp4a", "silva", "slc24a5", "snail2",
#"sox10", "sox9b", "tfap2a", "tfap2e", "tfec", "tyr", "tyrp1b")

Genes_dt <- Clust_dt

for (gene in genes) {
	gene_thres <- 700
	gene_n <- paste( "n_", gene, sep = "")
	gene_dt <- cells[X == gene]
	gene_dt[,X:=NULL]	
	Genes_dt[,(gene_n):=as.integer( t( gene_dt))]
	gene_f <- paste( "f_", gene, sep = "")
	Genes_dt[, (gene_f):= factor( as.integer( get(gene_n) > gene_thres ), levels = c(0,1)) ]
	gene_p <- paste( "p_", gene, sep = "")
	Gene_res<-Genes_dt[order(get(gene_f)), .( exp_i = levels( get( gene_f)), t = as.vector( table( get( gene_f)))/.N ), by = class][exp_i == 1][order(class)]
	Cl_stats[, (gene_p) := Gene_res[,t]]

	PlotFileName <- file.path( getwd(), "Plot", "Tanyas_genes", paste(Prefix, "_gene_", gene, "_distr.png" ,sep = ""))

	png(filename = PlotFileName)
		plot(Cl_stats[,X], Cl_stats[,Y], cex = 3, main = gene, xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min(Cl_stats[,Y]),max(Cl_stats[,Y])))
		points(Cl_stats[,X], Cl_stats[,Y], xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min( Cl_stats[,Y]),max(Cl_stats[,Y])), cex = 3 , 
				col = colors[1+99*round(Cl_stats[,get(gene_p)],digits = 2)], pch = 19)
		text(Cl_stats[,X], Cl_stats[,Y], xlim = c( min(Cl_stats[,X]),max(Cl_stats[,X])), ylim = c( min(Cl_stats[,Y]),max(Cl_stats[,Y])), labels = Cl_stats[,class])
		
	dev.off()	
	
			} #gene
 

