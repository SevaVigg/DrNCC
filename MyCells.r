require("Seurat")

#uses ipmc Seurat object created by seurat.r

dir.create(file.path( getwd(), "Plot", "VlnPlots_initial"), showWarnings = FALSE)
dir.create(file.path( getwd(), "Plot", "VlnPlots_ltk_P_mitfa_P"), showWarnings = FALSE)

MB <- ipmc@data[, FetchData(ipmc, "silva") < 2.8 & FetchData(ipmc, "mlphb") > 2]
ltk_pos_mitfa_pos <- ipmc@data[, FetchData(ipmc, "ltk") > 2.5 & FetchData(ipmc, "mitfa" > 2.5 ]

ipmc_my <- ipmc
ipmc_ltk_P_mitfa_P <- ipmc

ipmc_my@ident[FetchData(ipmc, "silva") < 2.8 & FetchData(ipmc, "mlphb") > 2] <- "MB"
ipmc_ltk_P_mitfa_P@ident[FetchData(ipmc,  "ltk") > 2.5 & FetchData(ipmc, "mitfa") > 2.5 ] <- "Ltk+Mitfa+"

genes <- rownames(ipmc_my@data)

for( gene in genes){
	VlnPlotFileName <- file.path(getwd(), "Plot", "VlnPlots_ltk_P_mitfa_P", paste(gene, ".png", sep =""))
	png(filename = VlnPlotFileName)
		VlnPlot(ipmc, gene)
	dev.off()
	}


