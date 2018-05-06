# This snippet calculates 3d tSNE object, which can be visualized by tSNE_3D snippet. Snippet uses _ipmc, the Seurat object 
# as a data source (see descriptions in seurat.r) 

require(Seurat)

ipmc <- ScaleData( object = ipmc, vars.to.regress = "nUMI") 
ipmc_nh <- RunTSNE( ipmc, dim_embed = 3, genes.use = rownames(ipmc@data) )



