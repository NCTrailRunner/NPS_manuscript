#srun --mem=100GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

set.seed(1)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(dplyr)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(glmGamPoi)
library(future)
library(harmony)
library(patchwork)

in.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/plots'
setwd(in.dir)
seu <- readRDS('6.integrate_datasets___method=harmony___ref=FALSE.umap.rds')
DefaultAssay(seu) <- "SCT"

# Split the object
seu_list <- SplitObject(object = seu, split.by = 'predicted_id')

# Loop through each split object
for (i in 1:length(seu_list)) {
  cell_type <- names(seu_list)[i]
  message(i, '/', length(seu_list), ': ', cell_type)
  
  # Get current object
  seu_i <- seu_list[[i]]
  
  # Save the original
  saveRDS(seu_i, paste0('7.split_object___', cell_type, '_.rds'))
  n <- ncol(seu_i[['harmony']]@cell.embeddings)
  n <- min(n, 30)  # Use at most 30 dimensions
  # Check if harmony reduction exists in this object
  if ("harmony" %in% names(seu_i@reductions)) {
    # Use harmony for clustering
    seu_i <- FindNeighbors(
	  seu_i, 
	  reduction = "harmony", 
	  dims = 1:n,  # Dynamically use the right number of dimensions
	  k.param = 50,  # Similar to k.nn in FindMultiModalNeighbors
	  prune.SNN = 0,  # Match the prune.SNN setting
	  graph.name = "wsnn"
    )
  } else {
    # Fallback to PCA if harmony is not available
    seu_i <- RunPCA(seu_i, assay = "SCT")
    seu_i <- FindNeighbors(
	  seu_i, 
	  reduction = "pca", 
	  dims = 1:n,  # Dynamically use the right number of dimensions
	  k.param = 50,  # Similar to k.nn in FindMultiModalNeighbors
	  prune.SNN = 0,  # Match the prune.SNN setting
	  graph.name = "wsnn"
    )
  }
  
  # Check if UMAP exists, if not run it
  if (!("umap" %in% names(seu_i@reductions))) {
    if ("harmony" %in% names(seu_i@reductions)) {
      seu_i <- RunUMAP(seu_i, reduction = "harmony", dims = 1:n)
    } else {
      seu_i <- RunUMAP(seu_i, reduction = "pca", dims = 1:n)
    }
  }
  
  # Generate plots at different resolutions
  setwd(out.dir)
  
  # Create plots at different resolutions
  resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7)
  
  for (res in resolutions) {
    res_str <- gsub("\\.", "", as.character(res))
    seu_i <- FindClusters(seu_i,graph.name = "wsnn", algorithm = 3, resolution = res)
    
    # Make sure the PNG device is properly closed after each plot
    filename <- paste0(cell_type, '.final.UMAP_r', res_str, '.png')
    png(filename, width = 700, height = 650)
    print(DimPlot(seu_i, reduction = "umap", label = TRUE, label.size = 10))
    dev.off()
  }
  
  # Save the final object with all clusters
  # setwd(in.dir)
  # saveRDS(seu_i, paste0('7.clusters___', cell_type, '.rds'))
}