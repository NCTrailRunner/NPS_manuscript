#srun --mem=100GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
#R

# Load necessary libraries
library(Seurat)
library(dplyr)

# Define paths
mayo_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/mayo.objects'  # path of the six objects
out_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects' # path of the integrated object

# Define the object paths
object_paths <- c(
  file.path(mayo_dir, "exc.clusters_resolution=0.3.mayo.annotated.rds"),
  file.path(mayo_dir, "inh.clusters_resolution=0.3.mayo.annotated.rds"),
  file.path(mayo_dir, "micro.clusters_resolution=0.1.mayo.annotated.rds"),
  file.path(mayo_dir, "astro.clusters_resolution=0.1.mayo.annotated.rds"),
  file.path(mayo_dir, "oligo.clusters_resolution=0.1.mayo.annotated.rds"),
  file.path(mayo_dir, "opc.clusters_resolution=0.1.mayo.annotated.rds")
)

# Load the integrated object
all <- readRDS(file.path(out_dir, '6.integrate_datasets___method=harmony___ref=FALSE.umap.rds'))

# Create a data frame to store the combined metadata
mayo_metadata <- data.frame(row.names = rownames(all@meta.data))

# Initialize columns with NA values
mayo_metadata$mayo.cell.type <- NA
mayo_metadata$mayo.cell.type.number <- NA

# Process each cell type object
for (obj_path in object_paths) {
  # Load the object
  cat("Processing:", obj_path, "\n")
  seurat_obj <- readRDS(obj_path)
  
  # Extract the metadata we need
  obj_metadata <- seurat_obj@meta.data[, c("mayo.cell.type", "mayo.cell.type.number"), drop = FALSE]
  
  # Get the cell IDs from this object
  cell_ids <- rownames(obj_metadata)
  
  # Update the mayo_metadata for cells present in this object
  matching_cells <- intersect(cell_ids, rownames(mayo_metadata))
  cat("Matching cells found:", length(matching_cells), "\n")
  
  if (length(matching_cells) > 0) {
    mayo_metadata[matching_cells, "mayo.cell.type"] <- obj_metadata[matching_cells, "mayo.cell.type"]
    mayo_metadata[matching_cells, "mayo.cell.type.number"] <- obj_metadata[matching_cells, "mayo.cell.type.number"]
  }
}

# Check how many cells were annotated
cat("Cells with mayo.cell.type annotation:", sum(!is.na(mayo_metadata$mayo.cell.type)), "out of", nrow(mayo_metadata), "\n")

# Add the mayo metadata to the integrated object
for (col in c("mayo.cell.type", "mayo.cell.type.number")) {
  all[[col]] <- mayo_metadata[[col]]
}

# Save the updated object
saveRDS(all, file.path(out_dir, '9.integrated_with_mayo_annotations.rds'))

# Print summary information
cat("Updated object saved with mayo annotations.\n")
cat("Dimensions of updated object:", dim(all), "\n")
cat("Metadata columns in updated object:", colnames(all@meta.data), "\n")