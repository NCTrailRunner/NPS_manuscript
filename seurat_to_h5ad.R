#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
# module load HDF5/1.14.3-rhel9
# seurat_to_h5ad.R - Convert Seurat object to h5ad format

# Command line arguments:
# 1. Input RDS file path
# 2. Output h5ad file path

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript seurat_to_h5ad.R <input_rds_file> <output_h5ad_file>")
}

input_file <- args[1]
output_file <- args[2]
# We want a standardized raw file name like "R6911631_raw.rds"
raw_file <- paste0(gsub("(_downsampled)?\\.h5ad$", "", output_file), "_raw.rds")

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(SeuratDisk)
library(dplyr)

# Function to rename cells:
# Remove the first token (assumed cell type) from the current name
# and then prepend the batch from orig.ident.
rename_cells <- function(seurat_obj) {
  old_names <- colnames(seurat_obj)
  batch <- as.character(seurat_obj@meta.data$batch)
  # Extract the last token from each cell name (assumed to be the barcode)
  barcode <- sapply(strsplit(old_names, "_"), function(x) tail(x, n = 1))
  new_names <- paste0(batch, "_", barcode)
  seurat_obj <- RenameCells(seurat_obj, new.names = new_names)
  return(seurat_obj)
}
cat("Loading Seurat object from", input_file, "\n")
seurat_obj <- readRDS(input_file)

if ("RNA" %in% names(seurat_obj@assays)) {
  DefaultAssay(seurat_obj) <- "RNA"
  cat("Set RNA as active assay\n")
  
  # Remove non-RNA assays.
  non_rna_assays <- setdiff(names(seurat_obj@assays), "RNA")
  for (assay_name in non_rna_assays) {
    seurat_obj[[assay_name]] <- NULL
  }
  cat("Removed non-RNA assays\n")
  
  # Force the "data" layer to equal the raw counts.
  seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]],
                                      layer = "data",
                                      new.data = seurat_obj[["RNA"]]@counts)
} else {
  stop("RNA assay not found in Seurat object")
}

cat("Cleaning cell names to include batch prefix...\n")
seurat_obj <- rename_cells(seurat_obj)

# Save the current "data" layer (raw counts) externally.
raw_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
saveRDS(raw_data, file = raw_file)
cat("Saved raw data layer separately to", raw_file, "\n")

# Now compute PCA so that the exported h5ad contains PCA embeddings.
if (!"pca" %in% names(seurat_obj@reductions)) {
  cat("Computing PCA for RNA assay\n")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = 30)
}

# Save cell count to a metadata file.
metadata_file <- gsub("\\.h5ad$", "_metadata.txt", output_file)
cat("Cell count:", ncol(seurat_obj), "\n")
write(ncol(seurat_obj), file = metadata_file)

# Save as h5Seurat (intermediate step).
h5seurat_file <- gsub("\\.h5ad$", ".h5Seurat", output_file)
cat("Saving as h5Seurat to", h5seurat_file, "\n")
SaveH5Seurat(seurat_obj, file = h5seurat_file, overwrite = TRUE)

cat("Converting to h5ad:", output_file, "\n")
Convert(h5seurat_file, dest = output_file, overwrite = TRUE)

if (file.exists(h5seurat_file)) {
  file.remove(h5seurat_file)
  cat("Removed temporary h5Seurat file\n")
}

rm(seurat_obj)
gc()

cat("Successfully converted to h5ad\n")
