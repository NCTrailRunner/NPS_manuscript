#srun --mem=100GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
#R

# --- Annotation Script ---
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

set.seed(1)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr) 
library(data.table)
library(tibble) 
library(scMayoMap)
library(patchwork)

in.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/mayo.objects'
#dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# Load the Mayo mapping
mayo.m1 <- read.csv('/hpc/group/adrc/zm77/redo0514/brain.mayo.m1.csv', header = TRUE)

# Get brain cell types from scMayoMap
brain_cell_types <- colnames(scMayoMapDatabase[scMayoMapDatabase$tissue == "brain", ])
brain.ct <- grep("brain:", brain_cell_types, value = TRUE)

# Define cell types and their optimal resolutions
cell_type_resolutions <- list(
  "astro" = 0.1,
  "exc" = 0.3,
  "inh" = 0.3,
  "micro" = 0.1,
  "oligo" = 0.1,
  "opc" = 0.1
)

# Function to add number suffixes to repetitive elements in a vector
add_suffix <- function(vec) {
  counts <- table(vec)
  suffixes <- ave(seq_along(vec), vec, FUN = function(x) seq(length(x)))
  result <- ifelse(counts[vec] > 1, paste0(vec, suffixes), paste0(vec, "1"))
  return(result)
}

# Process each cell type
for (ct in names(cell_type_resolutions)) {
  message(paste0("Processing cell type: ", ct))
  
  # Load the clustered Seurat object
  setwd(in.dir)
  seu_file <- paste0('7.clusters___', ct, '_resolution=', cell_type_resolutions[[ct]], '.rds')
  
  if (!file.exists(seu_file)) {
    message(paste0("File not found: ", seu_file))
    next
  }
  
  seurat.obj <- readRDS(seu_file)
  DefaultAssay(seurat.obj) <- "SCT"
  
  # Check how many clusters there are
  cluster_counts <- table(seurat.obj$seurat_clusters)
  num_clusters <- length(cluster_counts)
  message(paste0("Number of clusters detected: ", num_clusters))
  
  # Convert ct to match the format in mayo.m1
  if (ct == 'opc') {ct_mayo_format <- 'OPC'} else {
    ct_mayo_format <- toupper(substr(ct, 1, 1)) %>% paste0(tolower(substr(ct, 2, nchar(ct))))
  }
  # Get the Mayo cell types that correspond to this cell type (removing NAs)
  mayo_cell_types <- mayo.m1$brain.mayo[mayo.m1$brain.m1 == ct_mayo_format]
  mayo_cell_types <- mayo_cell_types[!is.na(mayo_cell_types)]
  
  if (length(mayo_cell_types) == 0) {
    message(paste0("No Mayo cell types mapped for ", ct))
    next
  }
  
  message(paste0("Found ", length(mayo_cell_types), " Mayo cell types for ", ct, ": ", paste(mayo_cell_types, collapse=", ")))
  
  # Check if there's only one Mayo cell type or multiple
  if (length(mayo_cell_types) == 1) {
    message(paste0("Using single Mayo cell type mapping for ", ct, ": ", mayo_cell_types))
    
    # For only one Mayo cell type
    metadata <- seurat.obj@meta.data
    clusters <- unique(metadata$seurat_clusters)
    cell.type <- strsplit(mayo_cell_types, ":")[[1]][2]
	if (cell.type == 'Oligodendrocyte precursor cell') {cell.type = 'OPC'}
    cell.type.df <- data.frame('seurat_clusters' = clusters, 'mayo.cell.type' = cell.type)
    
    # Add suffixes to distinguish multiple clusters with the same cell type
    mayo.ct <- cell.type.df$mayo.cell.type
    mayo.cluster <- add_suffix(mayo.ct)
    mayo.cluster <- as.character(mayo.cluster)
    cell.type.df$mayo.cell.type.number <- mayo.cluster
    
    # Add annotations to metadata
    metadata$row.names <- rownames(metadata)
    new.meta <- merge(metadata, cell.type.df, by = "seurat_clusters")
    rownames(new.meta) <- new.meta$row.names
    seurat.obj@meta.data <- new.meta
    
  } else {
    message(paste0("Using multiple Mayo cell types mapping for ", ct, ": ", paste(mayo_cell_types, collapse=", ")))
    
    # Create database filtered to relevant cell types
    brain.db <- scMayoMapDatabase[, c('tissue', 'gene', mayo_cell_types)]
    brain.db <- brain.db[brain.db$tissue == 'brain', ]
    
    # Different approach based on number of clusters
    if (num_clusters <= 1) {
      message("Only one cluster detected. Using top expressed genes approach.")
      
      # Get a representative signature for this cell type using top expressed genes
      avg_expr <- AverageExpression(seurat.obj, assays = "SCT", slot = "data")$SCT
      top_genes <- names(sort(avg_expr[,1], decreasing = TRUE))[1:200]  # Top 200 expressed genes
      
      # Create a pseudo-markers dataframe
      pseudo_markers <- data.frame(
        p_val = rep(0.001, length(top_genes)),
        avg_log2FC = seq(from = 5, to = 1, length.out = length(top_genes)),
        pct.1 = rep(0.9, length(top_genes)),
        pct.2 = rep(0.1, length(top_genes)),
        p_val_adj = rep(0.001, length(top_genes)),
        cluster = rep(0, length(top_genes)),
        gene = top_genes
      )
      
      # Use this for scMayoMap
      scMayoMap.obj <- scMayoMap(data = pseudo_markers, database = brain.db, tissue = 'brain')
      
    } else {
      # Standard approach for multiple clusters
      message("Finding markers for multiple clusters using MAST...")
      seurat.markers <- FindAllMarkers(seurat.obj, method = 'MAST')
      scMayoMap.obj <- scMayoMap(data = seurat.markers, database = brain.db, tissue = 'brain')
    }
    
    # Get cluster annotations
    mayo <- scMayoMap.obj
    clusters <- as.numeric(unique((mayo$markers$cluster)))
    
    # Use the cell type of the highest score for each cluster
    cell.type <- c()
    for (i in clusters) {
      cluster_celltypes <- mayo$markers$celltype[mayo$markers$cluster == i]
      if (length(cluster_celltypes) > 0) {
        cell.type <- c(cell.type, cluster_celltypes[length(cluster_celltypes)])
      } else {
        # Handle case where no cell type is assigned
        cell.type <- c(cell.type, paste0("Unknown_", i))
      }
    }
    
    # Create dataframe with cluster annotations
    cell.type.df <- data.frame('seurat_clusters' = clusters, 'mayo.cell.type' = cell.type)
    
    # Add suffixes for duplicate cell types
    mayo.ct <- cell.type.df$mayo.cell.type
    mayo.cluster <- add_suffix(mayo.ct)
    mayo.cluster <- as.character(mayo.cluster)
    cell.type.df$mayo.cell.type.number <- mayo.cluster
    
    # Add annotations to metadata
    metadata <- seurat.obj@meta.data
    metadata$row.names <- rownames(metadata)
    new.meta <- merge(metadata, cell.type.df, by = "seurat_clusters")
    rownames(new.meta) <- new.meta$row.names
    seurat.obj@meta.data <- new.meta
  }
  
  # Save the annotated object
  setwd(out.dir)
  saveRDS(seurat.obj, paste0(ct, '.clusters_resolution=', cell_type_resolutions[[ct]], '.mayo.annotated.rds'))
  
  # Create visualization of the annotations
  p1 <- DimPlot(seurat.obj, reduction = "umap", group.by = "seurat_clusters", 
                label = TRUE, repel = TRUE) + 
    ggtitle(paste0(ct, " - Clusters"))
  
  p2 <- DimPlot(seurat.obj, reduction = "umap", group.by = "mayo.cell.type", 
                label = TRUE, repel = TRUE) + 
    ggtitle(paste0(ct, " - Mayo Cell Types"))
  
  combined_plot <- p1 + p2
  
  # Save the visualization
  png(paste0(ct, ".mayo_annotation.png"), width = 1600, height = 800, res = 100)
  print(combined_plot)
  dev.off()
  
  message(paste0("Completed processing for ", ct))
}

message("All cell types processed. Annotated objects saved to: ", out.dir)