#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
# module load HDF5/1.14.3-rhel9
# process_single_sample.R
# This script processes a single sample from the h5Seurat files, using only RNA assay

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(SeuratDisk)
library(dplyr)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript process_single_sample.R <sample_id> <h5seurat_dir> <mapping_csv> <output_dir>")
}

sample_id <- args[1]
h5seurat_dir <- args[2]
mapping_csv <- args[3]
output_dir <- args[4]

cat("\nProcessing sample:", sample_id, "\n")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the mapping file
cat("Loading mapping file:", mapping_csv, "\n")
cell_data <- read.csv(mapping_csv)

# Standardize column names
if ("barcode" %in% colnames(cell_data)) {
  cell_id_col <- "barcode"
} else if ("Barcode" %in% colnames(cell_data)) {
  cell_id_col <- "Barcode"
  cell_data$barcode <- cell_data$Barcode
  cell_id_col <- "barcode"
}

if ("grouping.by" %in% colnames(cell_data)) {
  cell_type_col <- "grouping.by"
} else if ("cell.type" %in% colnames(cell_data)) {
  cell_type_col <- "cell.type"
}

if ("individualID" %in% colnames(cell_data)) {
  sample_id_col <- "individualID"
} else if ("IndividualID" %in% colnames(cell_data)) {
  sample_id_col <- "IndividualID"
  cell_data$individualID <- cell_data$IndividualID
  sample_id_col <- "individualID"
}

# Map cell types in the CSV to h5Seurat files
celltype_to_h5seurat <- list(
  "Microglia" = "microglia.h5Seurat",
  "Astrocyte" = "astrocytes.h5Seurat",
  "Oligodendrocytes" = "oligodendroglia.h5Seurat",
  "Inhibitory Neurons" = "inhibitory.h5Seurat",
  "Excitatory Neurons" = c("cux2+.h5Seurat", "cux2-.h5Seurat"),
  "OPCs" = "oligodendroglia.h5Seurat"  # Assuming OPCs are included in oligodendroglia
)

# Get cells for this sample
sample_cells <- cell_data[cell_data[[sample_id_col]] == sample_id, ]

if (nrow(sample_cells) == 0) {
  cat("  No cells found for sample", sample_id, "\n")
  quit(status = 0)
}

cat("  Found", nrow(sample_cells), "cells in mapping file\n")

# Create lists to store cells from each object
sample_data_list <- list()

# Process each cell type in this sample
for (cell_type in unique(sample_cells[[cell_type_col]])) {
  if (!cell_type %in% names(celltype_to_h5seurat)) {
    cat("  Unknown cell type:", cell_type, "- skipping\n")
    next
  }
  
  # Get cells for this sample and cell type
  cell_type_cells <- sample_cells[sample_cells[[cell_type_col]] == cell_type, ]
  
  if (nrow(cell_type_cells) == 0) {
    next
  }
  
  cat("  Processing", nrow(cell_type_cells), "cells from", cell_type, "\n")
  
  # Get h5Seurat file(s) for this cell type
  h5_files <- celltype_to_h5seurat[[cell_type]]
  
  # Handle multiple h5Seurat files for a cell type (e.g., Excitatory Neurons)
  cell_type_seurat_list <- list()
  
  for (h5_file in h5_files) {
    file_path <- file.path(h5seurat_dir, h5_file)
    
    if (!file.exists(file_path)) {
      cat("    File not found:", file_path, "\n")
      next
    }
    
    # Load h5Seurat file
    cat("    Loading", h5_file, "\n")
    tryCatch({
      seu <- LoadH5Seurat(file_path)
      
      # Set RNA as the active assay
      if ("RNA" %in% names(seu@assays)) {
        DefaultAssay(seu) <- "RNA"
        cat("    Set RNA as active assay\n")
      } else {
        cat("    Warning: RNA assay not found in", h5_file, "\n")
        next
      }
      
      # Find cells that exist in the Seurat object
      valid_cells <- intersect(cell_type_cells$barcode, colnames(seu))
      
      if (length(valid_cells) == 0) {
        cat("    No matching cells found in", h5_file, "\n")
        next
      }
      
      cat("    Found", length(valid_cells), "matching cells in", h5_file, "\n")
      
      # Subset the Seurat object
      seu_subset <- subset(seu, cells = valid_cells)
      
      # Add to list with a unique name based on h5_file
      h5_short_name <- sub("\\.h5Seurat$", "", h5_file)
      cell_type_seurat_list[[h5_short_name]] <- seu_subset
      
      # Clean up memory
      rm(seu)
      gc()
      
    }, error = function(e) {
      cat("    Error loading", file_path, ":", e$message, "\n")
    })
  }
  
  # Merge all objects for this cell type if needed
  if (length(cell_type_seurat_list) > 0) {
    if (length(cell_type_seurat_list) == 1) {
      sample_data_list[[cell_type]] <- cell_type_seurat_list[[1]]
    } else {
      # Ensure RNA is the active assay in all objects before merging
      for (i in 1:length(cell_type_seurat_list)) {
        if ("RNA" %in% names(cell_type_seurat_list[[i]]@assays)) {
          DefaultAssay(cell_type_seurat_list[[i]]) <- "RNA"
        }
      }
      
      # Merge multiple objects for this cell type (e.g., cux2+ and cux2-)
      cat("    Merging", length(cell_type_seurat_list), "objects for", cell_type, "\n")
      
      merged_cell_type <- merge(
        cell_type_seurat_list[[1]],
        y = cell_type_seurat_list[-1],
        add.cell.ids = names(cell_type_seurat_list)
      )
      sample_data_list[[cell_type]] <- merged_cell_type
    }
    
    # Add original cell type information
    sample_data_list[[cell_type]]$original_cell_type <- cell_type
  }
  
  # Clean up memory
  rm(cell_type_seurat_list)
  gc()
}

if (length(sample_data_list) == 0) {
  cat("  No valid cells found for sample", sample_id, "\n")
  quit(status = 0)
}

# Merge all cell types for this sample
cat("  Merging", length(sample_data_list), "cell types\n")

# Ensure the RNA assay is active in all objects before merging
for (i in 1:length(sample_data_list)) {
  if ("RNA" %in% names(sample_data_list[[i]]@assays)) {
    DefaultAssay(sample_data_list[[i]]) <- "RNA"
  }
}

if (length(sample_data_list) == 1) {
  merged_sample <- sample_data_list[[1]]
} else {
  merged_sample <- merge(
    sample_data_list[[1]], 
    y = sample_data_list[-1],
    add.cell.ids = names(sample_data_list)
  )
}

# Add sample ID to metadata
merged_sample$sample_id <- sample_id

# Add cell type information from original assignment
if (length(sample_data_list) > 1) {
  # Extract cell type from cell name prefix when multiple cell types
  cell_prefixes <- sapply(strsplit(colnames(merged_sample), "_"), function(x) x[1])
  merged_sample$cell_type <- cell_prefixes
} else {
  # Use the original cell type when only one cell type
  merged_sample$cell_type <- names(sample_data_list)[1]
}

# Make sure RNA is the active assay
if ("RNA" %in% names(merged_sample@assays)) {
  DefaultAssay(merged_sample) <- "RNA"
  
  # Check what's available in the RNA assay
  cat("  RNA assay contains:\n")
  for (layer_name in names(merged_sample@assays$RNA)) {
    cat("    ", layer_name, "\n")
  }
}

# Count cells
total_cells <- ncol(merged_sample)
cat("  Created merged object with", total_cells, "cells\n")

# Save merged sample
output_file <- file.path(output_dir, paste0(sample_id, ".rds"))
saveRDS(merged_sample, file = output_file)
cat("  Saved to", output_file, "\n")

# Also save a status file to track completed samples
status_file <- file.path(output_dir, "completed_samples.txt")
if (file.exists(status_file)) {
  completed_samples <- readLines(status_file)
} else {
  completed_samples <- character(0)
}
writeLines(c(completed_samples, sample_id), status_file)

cat("  Completed processing sample:", sample_id, "\n")