#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
#R

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

library(data.table)
library(ggplot2)
library('Seurat')
library('plyr')
library('dplyr')
library(ggvenn)

out.dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/summary/AD_DEP_ovl_Normal_DEP'

setwd ('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/nebula_results_AD_DEPR')
degs_list <- list.files(pattern = "_two_step_FDR.csv")

degs.all <- lapply(degs_list, 
			   function(x){
			   tmp <- read.csv(x)
			   tmp <- tmp[abs (tmp$log2FC) > 0.2, ]
			   tmp <- tmp[tmp$FDR_Adj_p_val < 0.1, ]
			   return(tmp)
			 })
names(degs.all) <- degs_list
degs.all <- degs.all[sapply(degs.all, nrow) > 0]

degs <- lapply (1:length (degs.all),
                function (x) {
                tmp <- degs.all[x]
				cluster <- gsub ('LOAD_DEPR_Normal_|_filtered_two_step_FDR.csv', '', names(tmp))
				tmp <- as.data.frame (degs.all[[x]])
				tmp$cell.type <- cluster
				rownames (tmp) <- NULL
				tmp <- tmp[order(tmp$FDR_Adj_p_val), ]
				return (tmp) 
				}) %>% bind_rows
degs.ad <- degs
degs.ad <- degs.ad %>%
  filter(is.finite(log2FC))
	
#degs.ad$cell.type <- recode (degs.ad$cell.type, 'Glutamatergic neuron+' = 'Glutamatergic neuron')
#celltype.RM <- unique (degs.ad$cell.type)

setwd (out.dir)
write.csv (degs.ad, 'AD.DEP.NEBULA.DEGs.two_step.csv')

setwd ('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/nebula_results_Normal_DEPR')
degs_list <- list.files(pattern = "_two_step_FDR.csv")

degs.all <- lapply(degs_list, 
			   function(x){
			   tmp <- read.csv(x)
			   tmp <- tmp[abs (tmp$log2FC) > 0.2, ]
			   tmp <- tmp[tmp$FDR_Adj_p_val < 0.1, ]
			   return(tmp)
			 })
names(degs.all) <- degs_list
degs.all <- degs.all[sapply(degs.all, nrow) > 0]

degs <- lapply (1:length (degs.all),
                function (x) {
                tmp <- degs.all[x]
				cluster <- gsub ('Normal_DEPR_|_filtered_two_step_FDR.csv', '', names(tmp))
				tmp <- as.data.frame (degs.all[[x]])
				tmp$cell.type <- cluster
				rownames (tmp) <- NULL
				tmp <- tmp[order(tmp$FDR_Adj_p_val), ]
				return (tmp) 
				}) %>% bind_rows
degs.mdd <- degs
degs.mdd <- degs.mdd %>%
  filter(is.finite(log2FC))
				
setwd (out.dir)
write.csv (degs.mdd, 'Normal_DEP.NEBULA.DEGs.two_step.csv')

ad.ct <- unique (degs.ad$cell.type)
mdd.ct <- unique (degs.mdd$cell.type)
cell.type <- intersect (ad.ct, mdd.ct)

# =============================================================================
# NEW SECTION: ALL DEGs (UP AND DOWN COMBINED) ANALYSIS
# =============================================================================
message("Analyzing all DEGs (up and down combined) between AD and Normal_DEP...")

inter.all <- list()
uniq.all <- list()
celltype.all <- list()

for (j in (1: length(cell.type)))
{
  # Get ALL DEGs (both up and down) for AD in this cell type
  set1.all <- degs.ad[degs.ad$cell.type == cell.type[j], ]$gene 
  
  # Get ALL DEGs (both up and down) for Normal_DEP in this cell type
  set2.all <- degs.mdd[degs.mdd$cell.type == cell.type[j], ]$gene 
  
  if (length(set1.all)==0 | length (set2.all)==0) {
    inter.all[[j]] <- 'NA'
  } else {
    inter.all[[j]] <- intersect(set1.all, set2.all)
  }
  
  print(paste("Cell type:", cell.type[j]))
  print(paste("AD_DEP DEGs:", length(set1.all)))
  print(paste("Normal_DEP DEGs:", length(set2.all)))
  print(paste("Common DEGs:", length(inter.all[[j]])))
  print("---")
  
  celltype.all[j] <- cell.type[j]
  
  # Calculate unique genes for each condition
  myl.all <- list(A = set1.all, B = set2.all)
  uniq.all[[j]] <- lapply(1:length(myl.all), function(n) setdiff(myl.all[[n]], unlist(myl.all[-n])))
  names(uniq.all[[j]]) <- list('AD_DEP', 'Normal_DEP')
  
  # Create Venn diagram for all DEGs
  setwd(out.dir)
  inter_list_all <- list('AD_DEP' = set1.all, 'Normal_DEP' = set2.all)
  venn.all <- ggvenn(inter_list_all, fill_color = c("purple", "blue"))
  
  # Save Venn diagram
  png(filename = paste0('AD_DEP_vs_Normal_DEP_', cell.type[j], '_', 'degs.venn.two_step.png'))
  print(venn.all)
  dev.off()
}

names(inter.all) <- celltype.all 
names(uniq.all) <- celltype.all

# Save common DEGs (all directions) to file
setwd(out.dir)
file.all <- "AD_DEP_vs_Normal_DEP.celltype.common.ALL.two_step.txt"
conn.all <- file(description = file.all, open = "w")

newlist.all <- lapply(seq_len(length(inter.all)), function(i){
  writeLines(text = paste0('intersect_', names(inter.all[i])), con = conn.all, sep = "\t")
  writeLines(text = paste(inter.all[[i]], collapse = ","), con = conn.all, sep = "\r") 
})

close(conn.all)

# Save unique DEGs (all directions) to file
output_file_all <- "LOAD_DEP_Normal_DEP.celltype.unique.ALL.two_step.txt"
file_conn_all <- file(output_file_all, open = "wt")

write_list_all <- function(data_list, file_connection) {
  for (top_level in names(data_list)) {
    writeLines(paste0("$", top_level), con = file_connection)
    sub_elements <- data_list[[top_level]]
    for (sub_category in names(sub_elements)) {
      writeLines(paste0("$", top_level, "$", sub_category), con = file_connection)
      genes <- sub_elements[[sub_category]]
      writeLines(paste(genes, collapse = ", "), con = file_connection)
    }
  }
}

write_list_all(uniq.all, file_conn_all)
close(file_conn_all)

# Create summary table for all DEGs analysis
summary_all <- data.frame(
  Cell_Type = names(inter.all),
  AD_DEGs = sapply(names(inter.all), function(ct) length(degs.ad[degs.ad$cell.type == ct, ]$gene)),
  Normal_DEP_DEGs = sapply(names(inter.all), function(ct) length(degs.mdd[degs.mdd$cell.type == ct, ]$gene)),
  Common_DEGs = sapply(inter.all, length),
  Unique_AD = sapply(uniq.all, function(x) length(x$AD)),
  Unique_Normal_DEP = sapply(uniq.all, function(x) length(x$Normal_DEP))
)

# Calculate percentages
summary_all$Percent_Common_of_AD <- round((summary_all$Common_DEGs / summary_all$AD_DEGs) * 100, 2)
summary_all$Percent_Common_of_Normal_DEP <- round((summary_all$Common_DEGs / summary_all$Normal_DEP_DEGs) * 100, 2)

# Save summary table
write.csv(summary_all, "LOAD.Normal_DEP.DEGs.summary.ALL.mayo.csv", row.names = FALSE)

message("Summary of ALL DEGs analysis:")
print(summary_all)

# =============================================================================
# ORIGINAL CODE: UP-REGULATED DEGs ANALYSIS
# =============================================================================

inter <- list()
uniq <- list()
celltype <- list()
uniq_name <- list()
for (j in (1: length(cell.type)))
{
  set1.up <- degs.ad[degs.ad$log2FC >0.2 & degs.ad$cell.type == cell.type[j] ,]$gene 
  
  
  set2.up <- degs.mdd[degs.mdd$log2FC >0.2 & degs.mdd$cell.type == cell.type[j] ,]$gene 
  
  
  if (length(set1.up)==0 | length (set2.up)==0) {
    inter[[j]] <- 'NA'} else {inter [[j]] <- intersect(set1.up,set2.up)}
  print (cell.type[j])
  print (length (inter[[j]]))
  celltype [j] <- cell.type[j]
  
  myl <- list(A = set1.up,
            B = set2.up
            )
  uniq[[j]] <- lapply(1:length(myl), function(n) setdiff(myl[[n]], unlist(myl[-n])))
  names (uniq[[j]]) <- list ('AD_DEP', 'Normal_DEP')
  setwd (out.dir)
  inter_list <-list('AD_DEP'=set1.up,'Normal_DEP'=set2.up)
  venn <- ggvenn(inter_list,  fill_color=c("purple", "blue"))
  png(filename = paste0('AD_DEP_vs_Normal_DEP_', cell.type[j], '_', 'degs.venn.celltype.up.two_step.png'))
  print(venn)
  #dev.copy(png,paste0('three_path_', cell.type[j], '_', 'degs.venn.png'))
  dev.off()
  
}


names(inter) <- celltype 
names (uniq) <- celltype

##deal with the celltype 'GABAergic neuron' where DEGs are more than 3000
# genes <- inter[[3]]
# setwd (ad.dir)
# data1 <- read.csv (paste0 ('AD.', cell.type[3], '_sig_genes_NEBULA.celltype.new.csv'))
# data1$ad.log2fc <- data1$log2FC
# data1 <- data1 [data1$gene %in% genes, c('gene', 'ad.log2fc')]

  
# setwd (mdd.dir)
# data2 <- read.csv (paste0 ('Normal_DEP.CONTROL.', cell.type[3], '_sig_genes_NEBULA.celltype.mayo.csv'))
# data2$mdd.log2fc <- data2$log2FC
# data2 <- data2 [data2$gene %in% genes, c('gene', 'mdd.log2fc')]

# inh <- full_join (data1, data2, by='gene')
# setwd (out.dir)
# #write.csv (inh, 'inh.up.common.csv')

# inh1 <- inh [inh$ad.log2fc >0.43 & inh$mdd.log2fc >0.43, ]
# nrow(inh1)
# inter[[3]] <- inh1$gene



setwd (out.dir)
file <- "AD_DEP_vs_Normal_DEP.celltype.common.up.two_step.txt"
conn <- file(description=file, open="w")

newlist <- lapply(seq_len(length(inter)), function(i){
    
     writeLines(text=paste0('intersect_', names(inter[i])), con=conn, sep="\t")
     writeLines(text=paste(inter[[i]], collapse=","), con=conn, sep="\r") 
   })

close(conn)

output_file <- "AD_DEP_vs_Normal_DEP.celltype.unique.up.two_step.txt"

# Open the file for writing
file_conn <- file(output_file, open = "wt")

# Function to write each top-level element and its sub-elements to the file
write_list <- function(data_list, file_connection) {
  # Iterate over the top-level elements (e.g., $Astrocyte, $Neuron)
  for (top_level in names(data_list)) {
    # Write the top-level element (e.g., $Astrocyte)
    writeLines(paste0("$", top_level), con = file_connection)
    
    # Get sub-elements (e.g., $Astrocyte$AD, $Astrocyte$Normal_DEP)
    sub_elements <- data_list[[top_level]]
    
    # Iterate over the sub-elements
    for (sub_category in names(sub_elements)) {
      # Write the sub-element name (e.g., $Astrocyte$AD)
      writeLines(paste0("$", top_level, "$", sub_category), con = file_connection)
      
      # Write the list of genes, comma-separated
      genes <- sub_elements[[sub_category]]
      writeLines(paste(genes, collapse = ", "), con = file_connection)
    }
  }
}

# Write the entire list to the file
write_list(uniq, file_conn)

# Close the file connection
close(file_conn)

# =============================================================================
# ORIGINAL CODE: DOWN-REGULATED DEGs ANALYSIS
# =============================================================================

inter <- list()
uniq <- list()
celltype <- list()
uniq_name <- list()
for (j in (1: length(cell.type)))
{
  set1.down <- degs.ad[degs.ad$log2FC < -0.2 & degs.ad$cell.type == cell.type[j] ,]$gene 
  
  
  set2.down <- degs.mdd[degs.mdd$log2FC < -0.2 & degs.mdd$cell.type == cell.type[j] ,]$gene #reverse the fold change direction based on our new understanding of NEBULA
  
  if (length(set1.down)==0 | length (set2.down)==0) {
    inter[[j]] <- 'NA'} else {inter [[j]] <- intersect(set1.down,set2.down)}
  print (cell.type[j])
  print (length (inter[[j]]))
  celltype [j] <- cell.type[j]
  
  myl <- list(A = set1.down,
            B = set2.down
           )
  uniq[[j]] <- lapply(1:length(myl), function(n) setdiff(myl[[n]], unlist(myl[-n])))
  names (uniq[[j]]) <- list ('AD_DEP', 'Normal_DEP')
  setwd (out.dir)
  inter_list <-list('AD_DEP'=set1.down,'Normal_DEP'=set2.down)
  venn <- ggvenn(inter_list,  fill_color=c("purple", "blue"))
  png(filename = paste0('AD_DEP_vs_Normal_DEP_', cell.type[j], '_', 'degs.venn.celltype.down.two_step.png'))
  print(venn)
  #dev.copy(png,paste0('three_path_', cell.type[j], '_', 'degs.venn.png'))
  dev.off()
  
}
names(inter) <- celltype 
names (uniq) <- celltype

setwd (out.dir)
file <- "AD_DEP_vs_Normal_DEP.celltype.common.down.two_step.txt"
conn <- file(description=file, open="w")

newlist <- lapply(seq_len(length(inter)), function(i){
      
     writeLines(text=paste0('intersect_', names(inter[i])), con=conn, sep="\t")
     writeLines(text=paste(inter[[i]], collapse=","), con=conn, sep="\r") 
   })

close(conn)

# Specify the output file path
output_file <- "AD_DEP_vs_Normal_DEP.celltype.unique.down.two_step.txt"

# Open the file for writing
file_conn <- file(output_file, open = "wt")

# Function to write each top-level element and its sub-elements to the file
write_list <- function(data_list, file_connection) {
  # Iterate over the top-level elements (e.g., $Astrocyte, $Neuron)
  for (top_level in names(data_list)) {
    # Write the top-level element (e.g., $Astrocyte)
    writeLines(paste0("$", top_level), con = file_connection)
    
    # Get sub-elements (e.g., $Astrocyte$AD, $Astrocyte$Normal_DEP)
    sub_elements <- data_list[[top_level]]
    
    # Iterate over the sub-elements
    for (sub_category in names(sub_elements)) {
      # Write the sub-element name (e.g., $Astrocyte$AD)
      writeLines(paste0("$", top_level, "$", sub_category), con = file_connection)
      
      # Write the list of genes, comma-separated
      genes <- sub_elements[[sub_category]]
      writeLines(paste(genes, collapse = ", "), con = file_connection)
    }
  }
}

# Write the entire list to the file
write_list(uniq, file_conn)

# Close the file connection
close(file_conn)

print(
  sessionInfo()
)