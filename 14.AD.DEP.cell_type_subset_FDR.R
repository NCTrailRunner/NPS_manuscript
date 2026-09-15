# Two-step FDR processing for all NEBULA results
# Step 1: Filter by |log2FC| > 0.2
# Step 2: Apply FDR correction to filtered genes

library(dplyr)

# Set working directory to your NEBULA results folder
setwd("/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/nebula_results_AD_DEPR")

# Function to process a single NEBULA file with two-step FDR
process_nebula_file <- function(file_path, log2fc_threshold = 0.2) {
  
  cat("Processing:", basename(file_path), "\n")
  
  # Read the original file
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Check if required columns exist
  required_cols <- c("gene", "log2FC", "p_val")
  if (!all(required_cols %in% colnames(data))) {
    cat("  Warning: Missing required columns in", basename(file_path), "\n")
    return(NULL)
  }
  
  # Step 1: Filter by absolute log2FC threshold
  data_filtered <- data[abs(data$log2FC) > log2fc_threshold, ]
  
  cat("  Original genes:", nrow(data), "\n")
  cat("  After |log2FC| >", log2fc_threshold, "filter:", nrow(data_filtered), "\n")
  
  # Step 2: Apply FDR correction to filtered subset
  if (nrow(data_filtered) > 0) {
    data_filtered$FDR_Adj_p_val <- p.adjust(data_filtered$p_val, method = "fdr")
    
    # Sort by FDR
    data_filtered <- data_filtered[order(data_filtered$FDR_Adj_p_val), ]
    
    # Count significant genes at different thresholds
    sig_005 <- sum(data_filtered$FDR_Adj_p_val < 0.05, na.rm = TRUE)
    sig_01 <- sum(data_filtered$FDR_Adj_p_val < 0.1, na.rm = TRUE)
    sig_02 <- sum(data_filtered$FDR_Adj_p_val < 0.2, na.rm = TRUE)
    
    cat("  Significant genes - FDR < 0.05:", sig_005, "\n")
    cat("  Significant genes - FDR < 0.1:", sig_01, "\n")
    cat("  Significant genes - FDR < 0.2:", sig_02, "\n")
    
    # Create output filename
    output_file <- gsub("_all_genes_NEBULA\\.csv$", "_filtered_two_step_FDR.csv", file_path)
    
    # Write the processed file
    write.csv(data_filtered, output_file, row.names = FALSE)
    cat("  Saved to:", basename(output_file), "\n\n")
    
    # Return summary statistics
    return(data.frame(
      File = basename(file_path),
      Original_genes = nrow(data),
      Filtered_genes = nrow(data_filtered),
      Sig_FDR_005 = sig_005,
      Sig_FDR_01 = sig_01,
      Sig_FDR_02 = sig_02,
      stringsAsFactors = FALSE
    ))
    
  } else {
    cat("  No genes passed the log2FC filter!\n\n")
    return(data.frame(
      File = basename(file_path),
      Original_genes = nrow(data),
      Filtered_genes = 0,
      Sig_FDR_005 = 0,
      Sig_FDR_01 = 0,
      Sig_FDR_02 = 0,
      stringsAsFactors = FALSE
    ))
  }
}

# Find all NEBULA result files
nebula_files <- list.files(pattern = "*all_genes_NEBULA\\.csv$", full.names = TRUE)

cat("Found", length(nebula_files), "NEBULA result files to process:\n")
for (i in seq_along(nebula_files)) {
  cat(i, ":", basename(nebula_files[i]), "\n")
}
cat("\n")

# Process all files
summary_list <- list()

for (i in seq_along(nebula_files)) {
  cat("=== Processing file", i, "of", length(nebula_files), "===\n")
  
  tryCatch({
    summary_stats <- process_nebula_file(nebula_files[i], log2fc_threshold = 0.2)
    if (!is.null(summary_stats)) {
      summary_list[[i]] <- summary_stats
    }
  }, error = function(e) {
    cat("Error processing", basename(nebula_files[i]), ":", e$message, "\n\n")
  })
}

# Combine all summary statistics
if (length(summary_list) > 0) {
  combined_summary <- do.call(rbind, summary_list)
  
  # Add cell type information
  combined_summary$Cell_Type <- gsub("LOAD_DEPR_Normal_|_all_genes_NEBULA\\.csv", "", combined_summary$File)
  
  # Reorder columns
  combined_summary <- combined_summary[, c("Cell_Type", "File", "Original_genes", "Filtered_genes", 
                                          "Sig_FDR_005", "Sig_FDR_01", "Sig_FDR_02")]
  
  # Save summary
  write.csv(combined_summary, "Two_Step_FDR_Processing_Summary.csv", row.names = FALSE)
  
  # Print summary
  cat("=== PROCESSING SUMMARY ===\n")
  print(combined_summary)
  
  # Create a more detailed summary
  cat("\n=== DETAILED SUMMARY ===\n")
  cat("Total files processed:", nrow(combined_summary), "\n")
  cat("Total original genes across all files:", sum(combined_summary$Original_genes), "\n")
  cat("Total filtered genes across all files:", sum(combined_summary$Filtered_genes), "\n")
  cat("Total significant genes (FDR < 0.05):", sum(combined_summary$Sig_FDR_005), "\n")
  cat("Total significant genes (FDR < 0.1):", sum(combined_summary$Sig_FDR_01), "\n")
  cat("Total significant genes (FDR < 0.2):", sum(combined_summary$Sig_FDR_02), "\n")
  
  # Calculate percentages
  pct_filtered <- round(sum(combined_summary$Filtered_genes) / sum(combined_summary$Original_genes) * 100, 2)
  cat("Percentage of genes passing |log2FC| > 0.2 filter:", pct_filtered, "%\n")
  
  if (sum(combined_summary$Filtered_genes) > 0) {
    pct_sig_005 <- round(sum(combined_summary$Sig_FDR_005) / sum(combined_summary$Filtered_genes) * 100, 2)
    pct_sig_01 <- round(sum(combined_summary$Sig_FDR_01) / sum(combined_summary$Filtered_genes) * 100, 2)
    cat("Percentage of filtered genes significant at FDR < 0.05:", pct_sig_005, "%\n")
    cat("Percentage of filtered genes significant at FDR < 0.1:", pct_sig_01, "%\n")
  }
  
} else {
  cat("No files were successfully processed!\n")
}

# Optional: Create visualization of results
if (exists("combined_summary") && nrow(combined_summary) > 0) {
  
  library(ggplot2)
  library(reshape2)
  
  # Prepare data for plotting
  plot_data <- combined_summary[, c("Cell_Type", "Sig_FDR_005", "Sig_FDR_01", "Sig_FDR_02")]
  plot_data_long <- melt(plot_data, id.vars = "Cell_Type", variable.name = "FDR_Threshold", value.name = "Significant_Genes")
  
  # Clean up FDR threshold labels
  plot_data_long$FDR_Threshold <- gsub("Sig_FDR_", "FDR < 0.", plot_data_long$FDR_Threshold)
  plot_data_long$FDR_Threshold <- factor(plot_data_long$FDR_Threshold, levels = c("FDR < 0.05", "FDR < 0.1", "FDR < 0.2"))
  
  # Create bar plot
  p <- ggplot(plot_data_long, aes(x = Cell_Type, y = Significant_Genes, fill = FDR_Threshold)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Two-Step FDR Results: Significant DEGs by Cell Type",
         subtitle = paste("Pre-filtered by |log2FC| > 0.2, then FDR correction applied"),
         x = "Cell Type", 
         y = "Number of Significant Genes",
         fill = "FDR Threshold") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(type = "seq", palette = "Blues")
  
  # Save plot
  ggsave("Two_Step_FDR_Results_Summary.png", p, width = 12, height = 8, dpi = 300)
  
  cat("\nSummary plot saved as: Two_Step_FDR_Results_Summary.png\n")
}

cat("\n=== TWO-STEP FDR PROCESSING COMPLETE ===\n")
cat("All processed files have '_filtered_two_step_FDR.csv' suffix\n")
cat("Summary saved as: Two_Step_FDR_Processing_Summary.csv\n")