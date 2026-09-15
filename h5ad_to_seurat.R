#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
# module load HDF5/1.14.3-rhel9
# h5ad_to_seurat.R - Convert h5ad file back to Seurat RDS without UMAP

# Command line arguments:
# 1. Input h5ad file path
# 2. Output RDS file path

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript h5ad_to_seurat.R <input_h5ad_file> <output_rds_file>")
}

input_file <- args[1]
output_file <- args[2]
# Always generate the raw file name as "RXXXXXX_raw.rds"
raw_file <- paste0(gsub("(_downsampled)?\\.h5ad$", "", input_file), "_raw.rds")

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

suppressPackageStartupMessages({
  library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
  library(SeuratDisk)
  library(dplyr)
  library(Matrix)
  library(reticulate)
})

# Force reticulate to use your conda environment.
python_path <- "/hpc/group/adrc/zm77/software/conda_envs/zmsc/bin/python"
reticulate::use_python(python_path, required = TRUE)
cat("Reticulate is using Python:", reticulate::py_config()$python, "\n")

# Function to rename cells:
# Remove the first token from the current cell name and prepend the batch (orig.ident).
rename_cells <- function(seurat_obj) {
  old_names <- colnames(seurat_obj)
  batch <- as.character(seurat_obj@meta.data$batch)
  # Extract the last token from each cell name (assumed to be the barcode)
  barcode <- sapply(strsplit(old_names, "_"), function(x) tail(x, n = 1))
  new_names <- paste0(batch, "_", barcode)
  seurat_obj <- RenameCells(seurat_obj, new.names = new_names)
  return(seurat_obj)
}

# Function to restore the raw slot from the externally saved raw counts.
# Instead of using a raw slot, we add a new assay "RNA_raw" that stores the raw counts.
restore_raw_assay <- function(seurat_obj, raw_file) {
  if (file.exists(raw_file)) {
    raw_data <- readRDS(raw_file)
    # Subset raw_data to include only the cells in the current object.
    common_cells <- intersect(colnames(raw_data), colnames(seurat_obj))
    if (length(common_cells) > 0) {
      raw_data <- raw_data[, common_cells, drop = FALSE]
      seurat_obj[["RNA_raw"]] <- CreateAssayObject(counts = raw_data)
      cat("Restored raw counts as assay 'RNA_raw' from", raw_file, "\n")
    } else {
      cat("No common cells found between current object and raw file\n")
    }
  } else {
    cat("Raw file", raw_file, "not found; raw assay not restored\n")
  }
  return(seurat_obj)
}


# --- Standard Conversion Function ---
standard_conversion <- function() {
  h5seurat_file <- gsub("\\.h5ad$", ".h5Seurat", input_file)
  cat("Attempting standard conversion:\nConverting to h5Seurat:", h5seurat_file, "\n")
  
  Convert(input_file, dest = "h5seurat", overwrite = TRUE, assay = "RNA")
  
  cat("Loading h5Seurat\n")
  seurat_obj <- LoadH5Seurat(h5seurat_file)
  
  if ("RNA" %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- "RNA"
    cat("Set RNA as active assay\n")
  } else {
    stop("RNA assay not found in the converted Seurat object")
  }
  
  seurat_obj <- rename_cells(seurat_obj)
  seurat_obj <- restore_raw_slot(seurat_obj, raw_file)
  
  cat("Saving Seurat object to", output_file, "\n")
  saveRDS(seurat_obj, file = output_file)
  
  if (file.exists(h5seurat_file)) {
    file.remove(h5seurat_file)
    cat("Removed temporary h5Seurat file\n")
  }
  return(TRUE)
}

# --- Alternative Conversion Function ---
alternative_conversion <- function() {
  cat("Attempting alternative conversion using reticulate and scanpy...\n")
  
  py_code <- "
import scanpy as sc
import pandas as pd

def read_h5ad(file_path):
    adata = sc.read_h5ad(file_path)
    if hasattr(adata, 'raw') and adata.raw is not None:
        raw = adata.raw
        if hasattr(raw.X, 'toarray'):
            counts = raw.X.toarray()
        else:
            counts = raw.X
        var_names = raw.var_names
    else:
        if hasattr(adata.X, 'toarray'):
            counts = adata.X.toarray()
        else:
            counts = adata.X
        var_names = adata.var_names
    df = pd.DataFrame(counts.T, index=var_names, columns=adata.obs_names)
    metadata = adata.obs.to_dict(orient='list')
    return {'df': df, 'metadata': metadata}
"
  cat("Running Python code to extract data from the h5ad file...\n")
  reticulate::py_run_string(py_code)
  
  read_h5ad <- reticulate::py$read_h5ad
  result <- read_h5ad(input_file)
  if (is.null(result)) {
    stop("Python function returned NULL")
  }
  
  counts_df <- result$df
  metadata <- result$metadata
  
  cat("Converting counts DataFrame to R matrix...\n")
  counts <- as.matrix(reticulate::py_to_r(counts_df))
  
  cat("Creating Seurat object using raw counts\n")
  seurat_obj <- CreateSeuratObject(counts = counts)
  
  if (!is.null(metadata)) {
    for (col in names(metadata)) {
      seurat_obj[[col]] <- metadata[[col]]
    }
  }
  
  seurat_obj <- rename_cells(seurat_obj)
  seurat_obj <- restore_raw_assay(seurat_obj, raw_file)
  
  cat("Saving Seurat object to", output_file, "\n")
  saveRDS(seurat_obj, file = output_file)
  return(TRUE)
}

conversion_success <- tryCatch({
  standard_conversion()
}, error = function(e) {
  cat("Standard conversion failed:", e$message, "\n")
  cat("Trying alternative conversion...\n")
  alt <- tryCatch({
    alternative_conversion()
  }, error = function(e2) {
    cat("Alternative conversion failed:", e2$message, "\n")
    FALSE
  })
  return(alt)
})

if (!conversion_success) {
  stop("All conversion methods failed")
} else {
  cat("Conversion completed successfully for", input_file, "\n")
}