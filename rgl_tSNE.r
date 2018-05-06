# This script requires MIGN object prepared by Seurat with three tSNE axes

require("Seurat")
require("rgl")
source("R/set_group_colors.r")

rgl.open()

x <- MIGN@dr$tsne@cell.embeddings[,1]
y <- MIGN@dr$tsne@cell.embeddings[,2]
z <- MIGN@dr$tsne@cell.embeddings[,3]


rgl.points(x,y,z, color = set_group_colors(MIGN@ident, group.col = c("black", "red", rep("grey", 11), "green", "grey")))





