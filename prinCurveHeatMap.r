if(!require(ComplexHeatmap)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("ComplexHeatmap")
}

tanyaGenes <- c("foxd3", "impdh1b",  "kita",     "ltk",      "mbpa", "mitfa", "mlphb", "oca2", "phox2b", "pnp4a", "silva", "slc24a5", "snail2",
                 "sox10", "sox9b", "tfec","tyrp1b",  "pax7b", "tfap2e", "tfap2a")

earlyId		<- colnames(LineageTree[[IP_linId]])[1]

CL3_linId	<- 5

IP_Id		<- colnames(LineageTree[[IP_linId]])[ length( colnames(LineageTree[[ IP_linId]])) ]
MC_Id		<- colnames(LineageTree[[MC_linId]])[ length( colnames(LineageTree[[ MC_linId]])) ]
CL3_Id		<- colnames(LineageTree[[CL3_linId]])[ length( colnames(LineageTree[[ CL3_linId]])) ]
MIXGN_Id	<- colnames(LineageTree[[IP_linId]])[ length( colnames(LineageTree[[ IP_linId]])) - 1 ]


curves_IP 	<- as.data.frame(logExps[tanyaGenes, names(curves[[ IP_linId]])])
curves_MC 	<- as.data.frame(logExps[tanyaGenes, names(curves[[ MC_linId]])])
curves_CL3 	<- as.data.frame(logExps[tanyaGenes, names(curves[[ CL3_linId]])])

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

CL3_clust 	<- ipmcMD@ident[names(curves_CL3)]
CL3_df 		<- data.frame( clust = CL3_clust)
CL3_cols		<- grey.colors( length(unique(CL3_df$clust)))

names(CL3_cols)	<- as.character(unique(CL3_df$clust))


both		<- intersect(names(MC_cols), names(IP_cols))
IP_cols[both]	<- terrain.colors( length(both))
MC_cols[both]	<- terrain.colors( length(both))
CL3_cols[both]	<- terrain.colors( length(both))

MC_cols[earlyId]	<- "red"
MC_cols[MC_Id]		<- "black"
MC_cols[MIXGN_Id]	<- "#FF99FF"
IP_cols[earlyId]	<- "red"
IP_cols[MIXGN_Id]	<- "#FF99FF"
IP_cols[IP_Id]		<- "cyan"
CL3_cols[earlyId]	<- "red"
CL3_cols[MIXGN_Id]	<- "#FF99FF"
CL3_cols[CL3_Id]	<- "orange"


png( file.path(heatmapPlotDir, "IP.png"), width = 800, height = 600)
plot(0,0, main = "Expression")

IP_topbar	<- columnAnnotation(IP_df, 
				    col = list( clust = IP_cols), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
					)

hplot <- Heatmap( 	curves_IP, 
			name = "log2Exp", 
			cluster_columns = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = IP_topbar,
			column_title = "IP_expression", 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
dev.off() 


png( file.path(heatmapPlotDir, "MC.png"), width = 800, height = 600)
plot(0,0, main = "Expression")

MC_topbar	<- columnAnnotation(MC_df, 
				    col = list( clust = MC_cols), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
					)

hplot <- Heatmap( 	curves_MC, 
			name = "log2Exp", 
			cluster_columns = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = MC_topbar,
			column_title = "MC_expression", 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
dev.off()

png( file.path(heatmapPlotDir, "CL3.png"), width = 800, height = 600)
plot(0,0, main = "Expression")

CL3_topbar	<- columnAnnotation(CL3_df, 
				    col = list( clust = CL3_cols), 
				    height = unit(30, "points"),
				    annotation_legend_param = list( nrow = 1)
					)

hplot <- Heatmap( 	curves_CL3, 
			name = "log2Exp", 
			cluster_columns = FALSE, 
			show_column_names = FALSE, 
			row_names_gp = gpar(fontsize = 24), 
			column_names_gp = gpar(fontsize = 10), 
			bottom_annotation = CL3_topbar,
			column_title = "CL3_expression", 
			clustering_distance_rows = "euclidean", 
			use_raster = TRUE, raster_device = "png", 
			)
draw(hplot, annotation_legend_side = "bottom")
dev.off() 

 
