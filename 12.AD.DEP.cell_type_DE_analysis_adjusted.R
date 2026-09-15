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
base_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out_dir <- file.path('/hpc/group/adrc/mwl17/nebula_results_AD_DEPR')
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
	
clinic <- clinic [clinic$ad == 'LOAD', ]
clinic$diagnosis <- paste0(clinic$ad, '_', clinic$depression)
clinic <- clinic [clinic$diagnosis == 'LOAD_DEPR' | clinic$diagnosis == 'LOAD_Normal', ]

print ('number of subjects with LOAD but without depression:')
print (nrow (clinic[clinic$diagnosis == 'LOAD_Normal', ]))
print ('number of subjects with LOAD and depression:')
print (nrow (clinic[clinic$diagnosis == 'LOAD_DEPR', ]))

# FIXED: Filter metadata first, then merge with clinical data
metadata <- metadata[metadata$DoubletFinder.score < 0.75, ]
metadata$rowname <- rownames(metadata)
metadata <- left_join(metadata, clinic, by='individualID')
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$rowname
metadata <- metadata [complete.cases(metadata$diagnosis), ]

# Remove cells with missing covariates
metadata <- metadata[!is.na(metadata$msex), ]
message("After removing cells with missing covariates: ", nrow(metadata), " cells")

# FIXED: Check if metadata rownames exist in Seurat object
cells_in_seurat <- colnames(seu)
cells_in_metadata <- rownames(metadata)
common_cells <- intersect(cells_in_seurat, cells_in_metadata)

message("Cells in Seurat object: ", length(cells_in_seurat))
message("Cells in filtered metadata: ", length(cells_in_metadata))
message("Common cells: ", length(common_cells))

if (length(common_cells) == 0) {
  stop("No common cells found between Seurat object and filtered metadata!")
}

# Keep only common cells
metadata <- metadata[common_cells, ]
message("Final metadata after matching with Seurat object: ", nrow(metadata), " cells")

# Create diagnosis variable (0 for LOAD with DEP, 1 for LOAD without DEP)
metadata$diagnosis <- ifelse(metadata$diagnosis == 'LOAD_DEPR', 1, 0)

# Update the Seurat object metadata (keep only filtered cells)
seu <- seu[, common_cells]  # FIXED: Use common_cells instead of rownames(metadata)
seu@meta.data <- metadata

# Get Mayo cell types and count how many cells in each
cell_types <- table(metadata$mayo.cell.type)
cell_types <- names(cell_types[cell_types >= 50])  # Only analyze cell types with at least 50 cells

message("Found ", length(cell_types), " Mayo cell types with at least 50 cells")
print(cell_types)

# Function to calculate adaptive FDR
adaptive_fdr <- function(p_values, alpha = 0.05) {
  m <- length(p_values)
  p_sorted <- sort(p_values)
  indices <- order(p_values)
  
  # Estimate proportion of true nulls (pi0)
  lambda <- seq(0, 0.95, 0.05)
  pi0_estimates <- sapply(lambda, function(l) {
    sum(p_values > l) / ((1 - l) * m)
  })
  pi0 <- min(pi0_estimates[pi0_estimates <= 1])
  pi0 <- max(pi0, 0.1)  # Don't let pi0 be too small
  
  # Calculate adaptive FDR
  fdr_vals <- numeric(m)
  for (i in 1:m) {
    fdr_vals[i] <- (pi0 * m * p_sorted[i]) / i
  }
  
  # Ensure monotonicity
  for (i in (m-1):1) {
    fdr_vals[i] <- min(fdr_vals[i], fdr_vals[i+1])
  }
  
  # Return in original order
  result_fdr <- numeric(m)
  result_fdr[indices] <- fdr_vals
  return(pmin(result_fdr, 1))
}

# Function to calculate two-stage FDR
two_stage_fdr <- function(p_values, alpha = 0.05) {
  # Stage 1: Apply FDR at alpha/2
  fdr_stage1 <- p.adjust(p_values, method = 'fdr')
  significant_stage1 <- which(fdr_stage1 <= alpha/2)
  
  if (length(significant_stage1) == 0) {
    return(fdr_stage1)
  }
  
  # Stage 2: Apply FDR only to significant genes from stage 1
  fdr_stage2 <- rep(1, length(p_values))
  fdr_stage2[significant_stage1] <- p.adjust(p_values[significant_stage1], method = 'fdr')
  
  return(fdr_stage2)
}

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
  
  # FIXED: Check that all cells exist in the Seurat object before subsetting
  cells_to_use <- rownames(meta)
  cells_available <- intersect(cells_to_use, colnames(seu))
  
  if (length(cells_available) != length(cells_to_use)) {
    message("Warning: Some cells not found in Seurat object. Using ", length(cells_available), " out of ", length(cells_to_use), " cells.")
    meta <- meta[cells_available, ]
  }
  
  # Get raw count data for these cells
  data.subset <- seu[["RNA_raw"]]$counts[, cells_available]
  
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
#  re <- nebula(count, sid, pred = df, offset = offsets, method = 'HL')
re <- nebula(count, sid, pred = df, offset = offsets, ncore = 1, cpc = 0, mincp = 0)

  
  # Process results
  result <- re$summary
  
  # Calculate log2FC (LOAD vs Normal)
  result$log2FC <- log2(exp(result$logFC_diagnosis))
  
  # ===========================================================================
  # MULTIPLE P-VALUE CORRECTION METHODS
  # ===========================================================================
  
  message("  - Applying multiple p-value correction methods...")
  
  # 1. Standard FDR (Benjamini-Hochberg)
  result$fdr_bh <- p.adjust(result$p_diagnosis, method = 'fdr', n = nrow(result))
  
  # 2. Bonferroni correction
  result$bonferroni <- p.adjust(result$p_diagnosis, method = 'bonferroni', n = nrow(result))
  
  # 3. Holm correction (step-down method)
  result$holm <- p.adjust(result$p_diagnosis, method = 'holm', n = nrow(result))
  
  # 4. Hochberg correction (step-up method) - RECOMMENDED
  result$hochberg <- p.adjust(result$p_diagnosis, method = 'hochberg', n = nrow(result))
  
  # 5. Hommel correction
  result$hommel <- p.adjust(result$p_diagnosis, method = 'hommel', n = nrow(result))
  
  # 6. Benjamini & Yekutieli (BY) correction
  result$by <- p.adjust(result$p_diagnosis, method = 'BY', n = nrow(result))
  
  # 7. Adaptive FDR
  result$adaptive_fdr <- adaptive_fdr(result$p_diagnosis)
  
  # 8. Two-stage FDR
  result$two_stage_fdr <- two_stage_fdr(result$p_diagnosis)
  
  # Select final columns for output
  result_final <- result[, c('gene', 'logFC_diagnosis', 'log2FC', 'p_diagnosis', 
                           'fdr_bh', 'bonferroni', 'holm', 'hochberg', 'hommel', 
                           'by', 'adaptive_fdr', 'two_stage_fdr')]
  
  # Rename columns for clarity
  colnames(result_final) <- c('gene', 'logFC', 'log2FC', 'p_val', 
                            'fdr_bh', 'bonferroni', 'holm', 'hochberg', 'hommel', 
                            'by', 'adaptive_fdr', 'two_stage_fdr')
  
  # Sort by original FDR
  result_final <- result_final[order(result_final$fdr_bh), ]
  
  # Count significant genes for each method
  sig_counts <- list(
    fdr_bh = sum(result_final$fdr_bh < 0.05, na.rm = TRUE),
    bonferroni = sum(result_final$bonferroni < 0.05, na.rm = TRUE),
    holm = sum(result_final$holm < 0.05, na.rm = TRUE),
    hochberg = sum(result_final$hochberg < 0.05, na.rm = TRUE),
    hommel = sum(result_final$hommel < 0.05, na.rm = TRUE),
    by = sum(result_final$by < 0.05, na.rm = TRUE),
    adaptive_fdr = sum(result_final$adaptive_fdr < 0.05, na.rm = TRUE),
    two_stage_fdr = sum(result_final$two_stage_fdr < 0.05, na.rm = TRUE)
  )
  
  message("  - Significant genes (p < 0.05) by method:")
  for (method in names(sig_counts)) {
    message(sprintf("    %s: %d genes", method, sig_counts[[method]]))
  }
  
  # Save comprehensive results
  setwd(out_dir)
  write.csv(result_final, paste0('LOAD_DEPR_Normal_', gsub(" ", "_", cell_type), '_allgenes_multiple_corrections.csv'), row.names = FALSE)
  
  # Save significant genes for each method (at p < 0.05)
  methods_to_save <- c('fdr_bh', 'hochberg', 'adaptive_fdr', 'two_stage_fdr')
  
  for (method in methods_to_save) {
    result_sig <- result_final[result_final[[method]] < 0.05, ]
    if (nrow(result_sig) > 0) {
      write.csv(result_sig, paste0('LOAD_DEPR_Normal_', gsub(" ", "_", cell_type), '_sig_genes_', method, '.csv'), row.names = FALSE)
    }
  }
  
  # Create comparison summary for this cell type
  comparison_df <- data.frame(
    Method = names(sig_counts),
    Significant_005 = unlist(sig_counts),
    Significant_01 = sapply(names(sig_counts), function(m) sum(result_final[[m]] < 0.1, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  write.csv(comparison_df, paste0('LOAD_DEPR_Normal_', gsub(" ", "_", cell_type), '_method_comparison.csv'), row.names = FALSE)
  
  # Create enhanced volcano plot with multiple methods
  if (nrow(result_final) > 10) {
    # Create volcano plot showing different significance thresholds
    result_plot <- result_final
    result_plot$significance <- "Not Significant"
    result_plot$significance[result_plot$fdr_bh < 0.05] <- "FDR < 0.05"
    result_plot$significance[result_plot$hochberg < 0.05] <- "Hochberg < 0.05"
    result_plot$significance[result_plot$adaptive_fdr < 0.05] <- "Adaptive FDR < 0.05"
    
    # Priority order for visualization
    result_plot$significance <- factor(result_plot$significance, 
                                     levels = c("Not Significant", "FDR < 0.05", "Hochberg < 0.05", "Adaptive FDR < 0.05"))
    
    p <- ggplot(result_plot, aes(x = log2FC, y = -log10(p_val))) +
      geom_point(aes(color = significance), alpha = 0.6, size = 1) +
      scale_color_manual(values = c("grey70", "blue", "orange", "red"), 
                        name = "Significance") +
      labs(title = paste0("Enhanced Volcano Plot: LOAD+DEP vs LOAD in ", cell_type),
           subtitle = paste0("Total genes: ", nrow(result_plot), 
                           " | FDR sig: ", sig_counts$fdr_bh,
                           " | Hochberg sig: ", sig_counts$hochberg,
                           " | Adaptive FDR sig: ", sig_counts$adaptive_fdr),
           x = "log2(Fold Change)",
           y = "-log10(p-value)") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5)
    
    ggsave(paste0('LOAD_DEPR_Normal_', gsub(" ", "_", cell_type), '_enhanced_volcano.png'), 
           p, width = 10, height = 8, dpi = 300)
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

# Create overall summary across all cell types
setwd(out_dir)
all_comparison_files <- list.files(pattern = "_method_comparison.csv")

if (length(all_comparison_files) > 0) {
  # Combine all comparison files
  all_comparisons <- lapply(all_comparison_files, function(f) {
    df <- read.csv(f)
    cell_type <- gsub("LOAD_DEPR_Normal_|_method_comparison.csv", "", f)
    df$Cell_Type <- cell_type
    return(df)
  })
  
  combined_comparison <- do.call(rbind, all_comparisons)
  write.csv(combined_comparison, "ALL_cell_types_method_comparison_summary.csv", row.names = FALSE)
  
  # Create summary plot
  library(ggplot2)
  library(reshape2)
  
  plot_data <- combined_comparison[, c("Cell_Type", "Method", "Significant_005")]
  plot_data_wide <- reshape(plot_data, idvar = "Cell_Type", timevar = "Method", direction = "wide")
  plot_data_melt <- melt(plot_data_wide, id.vars = "Cell_Type")
  plot_data_melt$Method <- gsub("Significant_005.", "", plot_data_melt$variable)
  
  p_summary <- ggplot(plot_data_melt, aes(x = Cell_Type, y = value, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Significant DEGs by Correction Method and Cell Type",
         x = "Cell Type", y = "Number of Significant Genes (p < 0.05)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(type = "qual", palette = "Set3")
  
  ggsave("ALL_cell_types_method_comparison_barplot.png", p_summary, width = 12, height = 8, dpi = 300)
}

# Summary of all analyses
summary_file <- file.path(out_dir, "nebula_analysis_AD_DEPR_enhanced_summary.txt")
analyzed_cell_types <- list.files(out_dir, pattern = "allgenes_multiple_corrections.csv")
sig_cell_types <- list.files(out_dir, pattern = "sig_genes_")

writeLines(c(
  paste("ENHANCED NEBULA Differential Expression Analysis Summary - AD+Depression Study"),
  paste("Date:", Sys.Date()),
  paste("Comparison: LOAD+Depression vs LOAD without Depression"),
  paste("Covariates: sex"),
  paste(""),
  paste("Multiple p-value correction methods applied:"),
  paste("- fdr_bh: Standard Benjamini-Hochberg FDR"),
  paste("- hochberg: Hochberg step-up method (RECOMMENDED)"),
  paste("- adaptive_fdr: Adaptive FDR (data-driven)"),
  paste("- two_stage_fdr: Two-stage FDR approach"),
  paste("- holm: Holm step-down method"),
  paste("- bonferroni: Conservative Bonferroni correction"),
  paste("- by: Benjamini-Yekutieli method"),
  paste("- hommel: Hommel closed testing procedure"),
  paste(""),
  paste("Total cell types analyzed:", length(analyzed_cell_types)),
  paste("Files generated per cell type:"),
  paste("- *_allgenes_multiple_corrections.csv: All genes with all correction methods"),
  paste("- *_sig_genes_[method].csv: Significant genes for each method"),
  paste("- *_method_comparison.csv: Summary statistics by method"),
  paste("- *_enhanced_volcano.png: Enhanced volcano plot"),
  paste(""),
  paste("Cell types analyzed:"),
  paste(" -", gsub("LOAD_DEPR_Normal_|\\_allgenes_multiple_corrections.csv", "", analyzed_cell_types))
), summary_file)

message("Enhanced analysis complete. Results saved to ", out_dir)
message("Check ALL_cell_types_method_comparison_summary.csv for overall comparison")

print(sessionInfo())