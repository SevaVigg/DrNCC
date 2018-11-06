if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}

tanyaGenes <- c("foxd3", "impdh1b",  "kita",     "ltk",      "mbpa", "mitfa", "mlphb", "oca2", "phox2b", "pnp4a", "silva", "slc24a5", "snail2",
                 "sox10", "sox9b", "tfec","tyrp1b",  "pax7b")

earlyId		<- colnames(LineageTree[[IP_linId]])[1]
IP_Id		<- colnames(LineageTree[[IP_linId]])[ length( colnames(LineageTree[[ IP_linId]])) ]
MC_Id		<- colnames(LineageTree[[MC_linId]])[ length( colnames(LineageTree[[ MC_linId]])) ]

curves_IP 	<- as.data.frame(logExps[tanyaGenes, names(curves[[ which( names( LineageTree) == IP_linId)]])])
curves_MC 	<- as.data.frame(logExps[tanyaGenes, names(curves[[ which( names( LineageTree) == MC_linId)]])])

greyCols	<- grey.colors( length( unique( ipmcMD@ident)))

heatmapPlotDir 	<- file.path(getwd(), "plot", "HeatMaps")
dir.create(heatmapPlotDir, showWarnings = FALSE)

IP_clust 	<- ipmcMD@ident[names(curves_IP)]
IP_df 		<- data.frame( clust = IP_clust)
IP_cols		<- grey.colors( length(unique(IP_df$clust)))

names(IP_cols)	<- as.character(unique(IP_df$clust))

MC_clust 	<- ipmcMD@ident[names(curves_MC)]
MC_df 		<- data.frame( clust = MC_clust)
MC_cols		<- grey.colors( length(unique(MC_df$clust)))

names(MC_cols)	<- as.character(unique(MC_df$clust))

both		<- intersect(names(MC_cols), names(IP_cols))
IP_cols[both]	<- terrain.colors( length(both))
MC_cols[both]	<- terrain.colors( length(both))

MC_cols[earlyId]	<- "red"
MC_cols[MC_Id]		<- "black"
IP_cols[earlyId]	<- "red"
IP_cols[IP_Id]		<- "cyan"


ha_topbar	<- columnAnnotation(IP_df, 
				    col = list( clust = IP_cols), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
					)


png( file.path(heatmapPlotDir, "IP.png"), width = 800, height = 600)
plot(0,0, main = "Expression")
hplot <- Heatmap( 	curves_IP, 
			name = "log2Exp", 
			cluster_columns = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = ha_topbar,
			column_title = "IP_expression", 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
dev.off() 
