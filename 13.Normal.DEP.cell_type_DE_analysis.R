#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/zmnebula
#R

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(data.table)
library(ggplot2)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(plyr)
library(dplyr)
library(stringr)
library(nebula)
library(missForest)
library(parallel)

# This can be run as an array job
i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()
# For testing, if not an array job, you can set i manually
if (is.na(i)) i <- 1

# Define directories
base_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out_dir <- file.path(base_dir, 'nebula_results_Normal_DEPR')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the integrated object with Mayo annotations
setwd(base_dir)
seu <- readRDS('9.integrated_with_mayo_annotations.rds')

# Make sure we have the raw RNA counts
if (!"RNA_raw" %in% names(seu@assays)) {
  stop("RNA_raw assay not found in the Seurat object. Please check that the raw counts are available.")
}

clinic <- read.csv("/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/ROSMAP_depr_clinical.csv")
# Get the metadata from Seurat object
metadata <- seu@meta.data
clinic <- clinic %>%
  mutate(depression = case_when(
    r_depres == 4 ~ "Normal",
    r_depres == 1 | r_depres == 2 | r_depres == 3 ~ "DEPR",
	TRUE ~ ""
	))

clinic <- clinic %>%
  mutate(ad = case_when(
    cogdx == 1 ~ "Normal",
    cogdx == 4 | cogdx == 5 ~ "LOAD",
	TRUE ~ ""
	))
	
clinic <- clinic [clinic$ad == 'Normal', ]
clinic$diagnosis <- paste0(clinic$ad, '_', clinic$depression)
clinic <- clinic [clinic$diagnosis == 'Normal_DEPR' | clinic$diagnosis == 'Normal_Normal', ]

print ('number of subjects with Normal but without depression:')
print (nrow (clinic[clinic$diagnosis == 'Normal_Normal', ]))
print ('number of subjects with Normal and depression:')
print (nrow (clinic[clinic$diagnosis == 'Normal_DEPR', ]))

metadata <- metadata[metadata$DoubletFinder.score < 0.75, ]
metadata$rowname <- rownames(metadata)
metadata <- left_join(metadata, clinic, by='individualID')
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$rowname
metadata <- metadata [complete.cases(metadata$diagnosis), ]

# Remove cells with missing covariates
metadata <- metadata[!is.na(metadata$msex), ]
message("After removing cells with missing covariates: ", nrow(metadata), " cells")

# Create diagnosis variable (0 for Normal without DEP(reference), 1 for Normal with DEP)
metadata$diagnosis <- ifelse(metadata$diagnosis == 'Normal_DEPR', 1, 0)

# Update the Seurat object metadata (keep only filtered cells)
seu <- seu[, rownames(metadata)]
seu@meta.data <- metadata

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
  
  # Get raw count data for these cells
  data.subset <- seu[["RNA_raw"]]$counts[, rownames(meta)]
  
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
  allgenes$sampID <- meta$individualID  # Using individualID for subject grouping
  allgenes$sex <- as.numeric(meta$msex)  # msex: 0 = female, 1 = male
  allgenes$agec <- meta$agec  # Age at death
  allgenes$pmi <- meta$pmi    # Post-mortem interval
  allgenes$nCount_RNA <- meta$nCount_RNA_raw
  allgenes$wellKey <- rownames(meta)
  rownames(allgenes) <- allgenes$wellKey
  
  # Reorder columns: metadata first, then genes
  metadata_cols <- c("diagnosis", "sampID", "sex", "agec", "pmi", "nCount_RNA", "wellKey")
  allgenes <- allgenes[, c(metadata_cols, gene_names)]
  
  # Handle missing values in covariates if any
  coldata <- allgenes[, 1:7]
  
  # Make sure sample_id is a factor
  coldata$sampID <- as.factor(coldata$sampID)
  
  # Filter for genes expressed in at least 10% of cells in either group
  message("  - Filtering genes by expression prevalence...")
  
  # Get cells for each condition
  load_cells <- as.matrix(t(allgenes[allgenes$diagnosis == 1, 8:ncol(allgenes)]))
  normal_cells <- as.matrix(t(allgenes[allgenes$diagnosis == 0, 8:ncol(allgenes)]))
  
  # Skip if one group has no cells
  if (ncol(load_cells) == 0 || ncol(normal_cells) == 0) {
    message("  - Not enough cells in one of the groups. Skipping.")
    return(NULL)
  }
  
  # Calculate percent of cells expressing each gene
  PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
  }
  
  pct.exp.load <- apply(X = load_cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  pct.exp.normal <- apply(X = normal_cells, MARGIN = 1, FUN = PercentAbove, threshold = 0)
  
  # Filter genes expressed in at least 10% of cells in either group
  alpha.min <- pmax(pct.exp.load, pct.exp.normal)
  genes.to.keep <- names(which(alpha.min >= 0.1))
  
  message("  - Keeping ", length(genes.to.keep), " genes expressed in at least 10% of cells")
  
  # Subset the data to these genes
  genedata <- data.subset[rownames(data.subset) %in% genes.to.keep, ]
  
  # Make sure gene data and cell data are in the same order
  genedata <- genedata[, rownames(coldata)]
  
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
  
  # Create design matrix with sex, age, and PMI as covariates
  df <- model.matrix(~ sex + diagnosis, data = coldata)
  
  message("  - Running NEBULA with ", nrow(count), " genes and ", ncol(count), " cells")
  
  # Run NEBULA
  re <- nebula(count, sid, pred = df, offset = offsets, method = 'HL')
  
  # Process results
  result <- re$summary
  
  # Calculate log2FC (LOAD vs Normal)
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
  write.csv(result, paste0('Normal_DEPR_', gsub(" ", "_", cell_type), '_all_genes_NEBULA.csv'), row.names = FALSE)
  
  if (nrow(result_sig) > 0) {
    write.csv(result_sig, paste0('Normal_DEPR_', gsub(" ", "_", cell_type), '_sig_genes_NEBULA.csv'), row.names = FALSE)
    message("  - Found ", nrow(result_sig), " differentially expressed genes at FDR < 0.05")
  } else {
    message("  - No significant differentially expressed genes found at FDR < 0.05")
  }
  
  # Create volcano plot
  if (nrow(result) > 10) {  # Only create plot if we have enough genes
    p <- ggplot(result, aes(x = log2FC, y = -log10(p_val))) +
      geom_point(aes(color = fdr < 0.05), alpha = 0.6) +
      scale_color_manual(values = c("grey", "red")) +
      labs(title = paste0("Volcano Plot: Normal_wt_DEP vs Normal_w_DEP ", cell_type),
           x = "log2(Fold Change)",
           y = "-log10(p-value)") +
      theme_minimal() +
      theme(legend.position = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed")
    
    ggsave(paste0('LOAD_DEPR_Normal_', gsub(" ", "_", cell_type), '_volcano.png'), p, width = 8, height = 6)
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
summary_file <- file.path(out_dir, "nebula_analysis_Normal_DEPR_summary.txt")
analyzed_cell_types <- list.files(out_dir, pattern = "all_genes_NEBULA.csv")
sig_cell_types <- list.files(out_dir, pattern = "sig_genes_NEBULA.csv")

writeLines(c(
  paste("NEBULA Differential Expression Analysis Summary - AD Study"),
  paste("Date:", Sys.Date()),
  paste("Comparison: Normal_wt_DEP vs Normal_w_DEP"),
  paste("Covariates: sex"),
  paste("Total cell types analyzed:", length(analyzed_cell_types)),
  paste("Cell types with significant genes:", length(sig_cell_types)),
  paste("Cell types analyzed:"),
  paste(" -", gsub("LOAD_vs_Normal_|\\_all_genes_NEBULA.csv", "", analyzed_cell_types))
), summary_file)

message("Analysis complete. Results saved to ", out_dir)

print(sessionInfo())