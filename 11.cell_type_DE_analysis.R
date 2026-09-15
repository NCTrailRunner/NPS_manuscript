#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/zmnebula
#R
source('/hpc/group/adrc/dcg27/african_american_multiome/scripts/config.R')
.libPaths(c('/hpc/group/adrc/dcg27/african_american_multiome/r_packages', .libPaths()))

set.seed(1)
library(nebula)
library(Matrix)
library(data.table)
library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(stringr)
library(missForest)
library(parallel)
library(future)

# This can be run as an array job
i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()
# For testing, if not an array job, you can set i manually
if (is.na(i)) i <- 1

# Define directories
base_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/MDD_new/objects'
out_dir <- file.path('/hpc/group/adrc/mwl17/nebula_results_yng')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the integrated object with Mayo annotations
setwd(base_dir)
seu <- readRDS('10.integrated_with_mayo_annotations.rds')

# Make sure we have the raw RNA counts
if (!"RNA_raw" %in% names(seu@assays)) {
  stop("RNA_raw assay not found in the Seurat object. Please check that the raw counts are available.")
}

# Get the metadata
metadata <- seu@meta.data

# Create diagnosis variable (1 for Case, 0 for Control)
metadata$diagnosis <- ifelse(metadata$Condition == 'Case', 1, 0)

# Get Mayo cell types and count how many cells in each
cell_types <- table(metadata$mayo.cell.type)
cell_types <- names(cell_types[cell_types >= 50])  # Only analyze cell types with at least 50 cells

message("Found ", length(cell_types), " Mayo cell types with at least 50 cells")
print(cell_types)

# Function to run NEBULA analysis for a given cell type
run_nebula_analysis <- function(cell_type) {
  message("Processing cell type: ", cell_type)
  
  # Subset metadata for this cell type
  meta <- metadata[metadata$mayo.cell.type == cell_type, ]
  
  # If not enough cells of this type, skip
  if (nrow(meta) < 50) {
    message("Not enough cells for cell type: ", cell_type, " (", nrow(meta), " cells). Skipping.")
    return(NULL)
  }
  
  # Check if we have both groups
  if (length(unique(meta$diagnosis)) < 2) {
    message("Only one diagnosis group found for cell type: ", cell_type, ". Skipping.")
    return(NULL)
  }
  
  # Get raw count data for these cells - only use cells that exist in both
  available_cells <- rownames(meta)[rownames(meta) %in% colnames(seu[["RNA_raw"]]$counts)]
  
  if (length(available_cells) < 50) {
    message("Not enough available cells for cell type: ", cell_type, " (", length(available_cells), " available). Skipping.")
    return(NULL)
  }
  
  # Update meta to only include available cells
  meta <- meta[available_cells, ]
  
  message("  - Using ", length(available_cells), " cells (", nrow(meta), " after filtering)")
  
  data.subset <- seu[["RNA_raw"]]$counts[, available_cells]
  
  # Filter genes with zero expression
  genes.use <- rowSums(data.subset) > 0
  genes.use <- names(genes.use[genes.use])
  
  message("  - Using ", length(genes.use), " genes with non-zero expression")
  
  # Subset the data to these genes
  data.subset <- data.subset[genes.use, ]
  
  # Prepare data for NEBULA
  allgenes <- t(as.data.frame(data.subset))
  allgenes <- as.data.frame(allgenes)
  
  # Add metadata
  ngenes <- ncol(allgenes)
  gene_names <- colnames(allgenes)  # Store original gene names
  
  allgenes$diagnosis <- meta$diagnosis
  allgenes$sampID <- meta$sample_id  # Using sample_id for subject grouping
  allgenes$sex <- as.numeric(meta$Sex == "Male")  # Convert sex to numeric (1 for Male, 0 for Female)
  allgenes$nCount_RNA <- meta$nCount_RNA_raw
  allgenes$wellKey <- rownames(meta)
  rownames(allgenes) <- allgenes$wellKey
  
  # Debug: Check the structure
  message("  - Number of metadata columns: ", 5)
  message("  - Number of gene columns: ", ngenes)
  message("  - Total columns in allgenes: ", ncol(allgenes))
  
  # Reorder columns: metadata first, then genes
  metadata_cols <- c("diagnosis", "sampID", "sex", "nCount_RNA", "wellKey")
  allgenes <- allgenes[, c(metadata_cols, gene_names)]
  
  # Handle missing values in covariates if any
  coldata <- allgenes[, 1:5]
  
  # Make sure sample_id is a factor
  coldata$sampID <- as.factor(coldata$sampID)
  
  # Filter for genes expressed in at least 10% of cells in either group
  message("  - Filtering genes by expression prevalence...")
  
  # Get cells for each condition - genes start from column 6
  gene_start_col <- 6
  case_cells <- as.matrix(t(allgenes[allgenes$diagnosis == 1, gene_start_col:ncol(allgenes)]))
  control_cells <- as.matrix(t(allgenes[allgenes$diagnosis == 0, gene_start_col:ncol(allgenes)]))
  
  # Skip if one group has no cells
  if (ncol(case_cells) == 0 || ncol(control_cells) == 0) {
    message("  - Not enough cells in one of the groups. Skipping.")
    return(NULL)
  }
  
  # Calculate percent of cells expressing each gene
  PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
  }
  
  pct.exp.case <- apply(X = case_cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  pct.exp.control <- apply(X = control_cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  
  # Filter genes expressed in at least 10% of cells in either group
  alpha.min <- pmax(pct.exp.case, pct.exp.control)
  genes.to.keep <- names(which(alpha.min >= 0.1))
  
  message("  - Keeping ", length(genes.to.keep), " genes expressed in at least 10% of cells")
  
  # Subset the data to these genes
  genedata <- data.subset[rownames(data.subset) %in% genes.to.keep, ]
  
  message("  - Gene data dimensions: ", nrow(genedata), " x ", ncol(genedata))
  message("  - Coldata dimensions: ", nrow(coldata), " x ", ncol(coldata))
  message("  - Column names match: ", all(colnames(genedata) %in% rownames(coldata)))
  
  # Make sure gene data and cell data are in the same order
  common_cells <- intersect(colnames(genedata), rownames(coldata))
  message("  - Common cells: ", length(common_cells))
  
  genedata <- genedata[, common_cells]
  coldata <- coldata[common_cells, ]
  
  # Order the data by sample ID
  coldata <- coldata[order(coldata$sampID), ]
  
  # Prepare for NEBULA
  count <- as.matrix(genedata)
  count <- count[, rownames(coldata)]
  
  # Set offset (library size normalization)
  offsets <- coldata$nCount_RNA
  offsets <- unname(offsets)
  
  # Create subject ID vector
  sid <- as.character(coldata$sampID)
  
  # Create design matrix with sex as covariate
  df <- model.matrix(~ sex + diagnosis, data = coldata)
  
  message("  - Running NEBULA with ", nrow(count), " genes and ", ncol(count), " cells")
  
  # Run NEBULA
  # re <- nebula(count, sid, pred = df, offset = offsets, method = 'HL')
  re <- nebula(count, sid, pred = df, offset = offsets, ncore = 1, cpc = 0, mincp = 0)
  

  # Process results
  result <- re$summary
  
  # Calculate log2FC (MDD vs Control)
  result$log2FC <- log2(exp(result$logFC_diagnosis))
  
  # Calculate FDR
  result$fdr <- p.adjust(result$p_diagnosis, method = 'fdr', n = nrow(result))
  
  # Rename columns
  result <- result[, c('gene', 'log2FC', 'p_diagnosis', 'fdr')]
  colnames(result) <- c('gene', 'log2FC', 'p_val', 'fdr')
  
  # Sort by FDR
  result <- result[order(result$fdr), ]
  
  # Extract significant genes
  result_sig <- result[which(result$fdr < 0.05), ]
  
  # Save results
  setwd(out_dir)
  write.csv(result, paste0('MDD_vs_Control_', gsub(" ", "_", cell_type), '_all_genes_NEBULA.csv'), row.names = FALSE)
  
  if (nrow(result_sig) > 0) {
    write.csv(result_sig, paste0('MDD_vs_Control_', gsub(" ", "_", cell_type), '_sig_genes_NEBULA.csv'), row.names = FALSE)
    message("  - Found ", nrow(result_sig), " differentially expressed genes at FDR < 0.05")
  } else {
    message("  - No significant differentially expressed genes found at FDR < 0.05")
  }
  
  # Create volcano plot
  if (nrow(result) > 10) {  # Only create plot if we have enough genes
    p <- ggplot(result, aes(x = log2FC, y = -log10(p_val))) +
      geom_point(aes(color = fdr < 0.05), alpha = 0.6) +
      scale_color_manual(values = c("grey", "red")) +
      labs(title = paste0("Volcano Plot: MDD vs Control in ", cell_type),
           x = "log2(Fold Change)",
           y = "-log10(p-value)") +
      theme_minimal() +
      theme(legend.position = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed")
    
    ggsave(paste0('MDD_vs_Control_', gsub(" ", "_", cell_type), '_volcano.png'), p, width = 8, height = 6)
  }
  
  # Return success
  return(TRUE)
}

# If running as array job, process the i-th cell type
if (!is.na(i) && i <= length(cell_types)) {
  result <- run_nebula_analysis(cell_types[i])
  if (is.null(result)) {
    message("Analysis failed or was skipped for cell type: ", cell_types[i])
  }
} else {
  # If not an array job or for testing, process all cell types sequentially
  message("Processing all cell types sequentially...")
  for (cell_type in cell_types) {
    result <- run_nebula_analysis(cell_type)
    gc()  # Clean up memory between runs
  }
}

# Summary of all analyses
summary_file <- file.path(out_dir, "nebula_analysis_summary.txt")
analyzed_cell_types <- list.files(out_dir, pattern = "all_genes_NEBULA.csv")
sig_cell_types <- list.files(out_dir, pattern = "sig_genes_NEBULA.csv")

writeLines(c(
  paste("NEBULA Differential Expression Analysis Summary"),
  paste("Date:", Sys.Date()),
  paste("Total cell types analyzed:", length(analyzed_cell_types)),
  paste("Cell types with significant genes:", length(sig_cell_types)),
  paste("Cell types analyzed:"),
  paste(" -", gsub("MDD_vs_Control_|\\_all_genes_NEBULA.csv", "", analyzed_cell_types))
), summary_file)

message("Analysis complete. Results saved to ", out_dir)

print(sessionInfo())