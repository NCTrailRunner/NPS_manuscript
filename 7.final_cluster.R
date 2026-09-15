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
library(ggplot2) 

in.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/plots'

# Define optimal resolutions for each cell type based on your evaluation
cell_type_resolutions <- list(
  "astro" = 0.1,
  "exc" = 0.3,
  "inh" = 0.3,
  "micro" = 0.1,
  "oligo" = 0.1,
  "opc" = 0.1
)

setwd(in.dir)
seu <- readRDS('6.integrate_datasets___method=harmony___ref=FALSE.umap.rds')
DefaultAssay(seu) <- "SCT"

# Split the object
seu_list <- SplitObject(object = seu, split.by = 'predicted_id')

# Loop through each split object
for (i in 1:length(seu_list)) {
  cell_type <- names(seu_list)[i]
  message(i, '/', length(seu_list), ': ', cell_type)
  
  # Skip if not one of our target cell types
  if (!cell_type %in% names(cell_type_resolutions)) {
    message("Skipping cell type: ", cell_type, " (not in target list)")
    next
  }
  
  # Get current object
  seu_i <- seu_list[[i]]
  setwd(in.dir)
  # Save the original split object
  saveRDS(seu_i, paste0('7.split_object___', cell_type, '_.rds'))
  
  # Get dimensions for harmony
  n <- ncol(seu_i[['harmony']]@cell.embeddings)
  n <- min(n, 30)  # Use at most 30 dimensions
  
  # Check if harmony reduction exists in this object
  if ("harmony" %in% names(seu_i@reductions)) {
    # Use harmony for clustering
    seu_i <- FindNeighbors(
      seu_i, 
      reduction = "harmony", 
      dims = 1:n,  
      k.param = 50,  
      prune.SNN = 0,  
      graph.name = "wsnn"
    )
  } else {
    # Fallback to PCA if harmony is not available
    seu_i <- RunPCA(seu_i, assay = "SCT")
    seu_i <- FindNeighbors(
      seu_i, 
      reduction = "pca", 
      dims = 1:n,  
      k.param = 50,  
      prune.SNN = 0,  
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
  
  # Get the optimal resolution for this cell type
  resolution <- cell_type_resolutions[[cell_type]]
  
  # Run clustering with the optimal resolution
  seu_i <- FindClusters(seu_i, graph.name = "wsnn", algorithm = 3, resolution = resolution)
  
  # Generate final plot with the optimal resolution
  setwd(out.dir)
  res_str <- gsub("\\.", "", as.character(resolution))
  filename <- paste0(cell_type, '.final.UMAP_r', res_str, '.optimal.png')
  png(filename, width = 700, height = 650)
  
  # Create and print the plot correctly
  p <- DimPlot(seu_i, reduction = "umap", label = TRUE, label.size = 10) + 
       ggtitle(paste0(cell_type, " - Resolution: ", resolution))
  print(p)
  
  dev.off()
  
  # Also create a UMAP with cells colored by sample_id to check batch effects
  if ("sample_id" %in% colnames(seu_i@meta.data)) {
    filename <- paste0(cell_type, '.sample_id.UMAP.png')
    png(filename, width = 700, height = 650)
    p_sample <- DimPlot(seu_i, reduction = "umap", group.by = "sample_id") + 
                ggtitle(paste0(cell_type, " - Colored by sample_id"))
    print(p_sample)
    dev.off()
  }
  
  # Save the final object with optimal clustering
  setwd(in.dir)
  saveRDS(seu_i, paste0('7.clusters___', cell_type, '_resolution=', resolution, '.rds'))
  
  # Also save just the cluster assignments for later reference
  clusters <- seu_i$seurat_clusters
  saveRDS(clusters, paste0('7.cluster_assignments___', cell_type, '_resolution=', resolution, '.rds'))
}

# Create a combined overview plot of all cell types
setwd(in.dir)
plot_list <- list()

for (cell_type in names(cell_type_resolutions)) {
  resolution <- cell_type_resolutions[[cell_type]]
  seu_file <- paste0('7.clusters___', cell_type, '_resolution=', resolution, '.rds')
  
  if (file.exists(seu_file)) {
    seu_i <- readRDS(seu_file)
    p <- DimPlot(seu_i, reduction = "umap", label = TRUE, 
               label.size = 6, repel = TRUE) + 
         ggtitle(paste0(cell_type, " (res:", resolution, ")")) +
         theme(plot.title = element_text(size = 12))
    plot_list[[cell_type]] <- p
  }
}

# Combine all plots if we have any
if (length(plot_list) > 0) {
  setwd(out.dir)
  # Calculate grid dimensions
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  # Combine plots
  combined_plot <- wrap_plots(plot_list, ncol = n_cols)
  
  # Save the combined plot
  png("all_cell_types_optimal_clusters.png", width = n_cols * 800, height = n_rows * 700, res = 100)
  print(combined_plot)
  dev.off()
}

message("Processing complete. Optimal clustering saved for each cell type at their respective resolutions:")
for (cell_type in names(cell_type_resolutions)) {
  message(paste0("- ", cell_type, ": resolution = ", cell_type_resolutions[[cell_type]]))
}