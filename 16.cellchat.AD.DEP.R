# srun --mem=200GB --pty bash -i
# conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
# R

source('/hpc/group/adrc/dcg27/african_american_multiome/scripts/config.R')
.libPaths(c('/hpc/group/adrc/dcg27/african_american_multiome/r_packages', .libPaths()))
set.seed(1)
library(Matrix)
library(data.table)
library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
options(stringsAsFactors = FALSE)
library (CellChat)
library(NMF)
library(ggalluvial)

# Define directories
base_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects'
out_dir <- file.path('/hpc/group/adrc/mwl17/cellchat_AD_DEPR')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the integrated object with Mayo annotations
setwd(base_dir)
combined.integrated <- readRDS('9.integrated_with_mayo_annotations.rds')

# Set default assay to SCT (which contains normalized data)
DefaultAssay(combined.integrated) <- "SCT"

# Load clinical data and process
clinic <- read.csv("/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/ROSMAP_depr_clinical.csv")

# Create depression and AD variables
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

# Filter for LOAD subjects only
clinic <- clinic[clinic$ad == 'LOAD', ]
clinic$diagnosis <- paste0(clinic$ad, '_', clinic$depression)
clinic <- clinic[clinic$diagnosis == 'LOAD_DEPR' | clinic$diagnosis == 'LOAD_Normal', ]

print('Clinical data diagnosis distribution:')
print(table(clinic$diagnosis))
print('Number of subjects with LOAD but without depression:')
print(nrow(clinic[clinic$diagnosis == 'LOAD_Normal', ]))
print('Number of subjects with LOAD and depression:')
print(nrow(clinic[clinic$diagnosis == 'LOAD_DEPR', ]))

# Get metadata from Seurat object
metadata <- combined.integrated@meta.data

# Filter by DoubletFinder score
metadata <- metadata[metadata$DoubletFinder.score < 0.75, ]
message("After DoubletFinder filtering: ", nrow(metadata), " cells")

# Save rownames before merge
metadata$cell_id <- rownames(metadata)

# Merge with clinical data
metadata <- left_join(metadata, clinic[, c('individualID', 'diagnosis', 'ad', 'depression', 'msex')], by='individualID')

# Restore rownames
rownames(metadata) <- metadata$cell_id

# Filter for cells that have clinical data (diagnosis not NA)
metadata_with_clinic <- metadata[!is.na(metadata$diagnosis), ]
message("Cells with clinical data: ", nrow(metadata_with_clinic), " cells")

# Remove cells with missing covariates
metadata_with_clinic <- metadata_with_clinic[!is.na(metadata_with_clinic$msex), ]
message("After removing cells with missing covariates: ", nrow(metadata_with_clinic), " cells")

# Check if we have cells for both conditions
print("Diagnosis distribution in cells with clinical data:")
print(table(metadata_with_clinic$diagnosis))

if (nrow(metadata_with_clinic) == 0) {
  stop("No cells found with matching clinical data!")
}

# Check cell type distribution by diagnosis
print("Cell type distribution by diagnosis:")
print(table(metadata_with_clinic$mayo.cell.type, metadata_with_clinic$diagnosis))

# Subset Seurat object to only cells with clinical data
cells_to_keep <- rownames(metadata_with_clinic)
combined.integrated <- combined.integrated[, cells_to_keep]

# Update metadata in Seurat object
combined.integrated@meta.data <- metadata_with_clinic

message("Final dataset: ", ncol(combined.integrated), " cells")

# Get normalized data from SCT assay for CellChat
data.input <- combined.integrated@assays$SCT@data

# Verify data structure
print("Data input dimensions:")
print(dim(data.input))
print("Data input class:")
print(class(data.input))

# Extract cells for each condition
meta <- combined.integrated@meta.data

cell.use1 <- rownames(meta)[meta$diagnosis == "LOAD_DEPR"] # LOAD with depression
meta1 <- meta[cell.use1, ]
print(paste("LOAD_DEPR cells:", length(cell.use1)))
print("LOAD_DEPR cell types:")
print(table(meta1$mayo.cell.type))

cell.use2 <- rownames(meta)[meta$diagnosis == "LOAD_Normal"] # LOAD without depression  
meta2 <- meta[cell.use2, ]
print(paste("LOAD_Normal cells:", length(cell.use2)))
print("LOAD_Normal cell types:")
print(table(meta2$mayo.cell.type))

# Find common cell types between both groups (excluding NA)
cell.type1 <- unique(meta1$mayo.cell.type)
cell.type1 <- cell.type1[!is.na(cell.type1)]
cell.type2 <- unique(meta2$mayo.cell.type) 
cell.type2 <- cell.type2[!is.na(cell.type2)]

cell.type <- intersect(cell.type1, cell.type2)
print("Common cell types between LOAD_DEPR and LOAD_Normal:")
print(cell.type)

# Check cell counts for each cell type in both conditions
if(length(cell.type) > 0) {
  cell_type_counts <- data.frame(
    CellType = cell.type,
    LOAD_DEPR = sapply(cell.type, function(x) sum(meta1$mayo.cell.type == x, na.rm = TRUE)),
    LOAD_Normal = sapply(cell.type, function(x) sum(meta2$mayo.cell.type == x, na.rm = TRUE))
  )
  print("Cell type counts by condition:")
  print(cell_type_counts)
  
  # Filter cell types that have at least 10 cells in both conditions (CellChat minimum)
  sufficient_cell_types <- cell_type_counts$CellType[
    cell_type_counts$LOAD_DEPR >= 10 & cell_type_counts$LOAD_Normal >= 10
  ]
  print(paste("Cell types with sufficient cells (>=10 in both conditions):", length(sufficient_cell_types)))
  print(sufficient_cell_types)
  
  # Update cell.type to only include those with sufficient cells
  cell.type <- sufficient_cell_types
  
  if(length(cell.type) == 0) {
    stop("No cell types have sufficient cells (>=10) in both conditions for CellChat analysis!")
  }
} else {
  stop("No common cell types found between conditions!")
}

###################################################################
# Analysis for LOAD_DEPR (LOAD with depression)
###################################################################

# Extract cells for LOAD_DEPR condition with sufficient cell types
cell.use1 <- rownames(meta)[meta$diagnosis == "LOAD_DEPR" & meta$mayo.cell.type %in% cell.type]
data.input1 <- data.input[, cell.use1]
metadata1 <- data.frame(labels = meta[cell.use1, "mayo.cell.type"], row.names = cell.use1)

print(paste("LOAD_DEPR cells for CellChat:", length(cell.use1)))
print("LOAD_DEPR cell type distribution for CellChat:")
print(table(metadata1$labels))

# Create CellChat object for LOAD_DEPR
cellchat1 <- createCellChat(object = data.input1, meta = metadata1, group.by = "labels")
cellchat1 <- addMeta(cellchat1, meta = metadata1)
cellchat1 <- setIdent(cellchat1, ident.use = "labels")
levels(cellchat1@idents) # show factor levels of the cell labels
groupSize1 <- as.numeric(table(cellchat1@idents)) # number of cells in each cell group

# Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
cellchat1@DB <- CellChatDB

# Preprocessing
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
# cellchat1 <- projectData(cellchat1, PPI.human)
cellchat1 <- smoothData(cellchat1, adj=PPI.human)

# Compute communication probabilities
cellchat1 <- computeCommunProb(cellchat1, raw.use = FALSE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
df.net1 <- subsetCommunication(cellchat1)

# Compute communication probabilities at pathway level
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat1 <- aggregateNet(cellchat1)
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")

# Save results
setwd(out_dir)
saveRDS(cellchat1, 'LOAD_DEPR.cellchat.rds')

###################################################################
# Analysis for LOAD_Normal (LOAD without depression)
###################################################################

# Extract cells for LOAD_Normal condition with sufficient cell types
cell.use2 <- rownames(meta)[meta$diagnosis == "LOAD_Normal" & meta$mayo.cell.type %in% cell.type]
data.input2 <- data.input[, cell.use2]
metadata2 <- data.frame(labels = meta[cell.use2, "mayo.cell.type"], row.names = cell.use2)

print(paste("LOAD_Normal cells for CellChat:", length(cell.use2)))
print("LOAD_Normal cell type distribution for CellChat:")
print(table(metadata2$labels))

# Create CellChat object for LOAD_Normal
cellchat2 <- createCellChat(object = data.input2, meta = metadata2, group.by = "labels")
cellchat2 <- addMeta(cellchat2, meta = metadata2)
cellchat2 <- setIdent(cellchat2, ident.use = "labels")
levels(cellchat2@idents) # show factor levels of the cell labels
groupSize2 <- as.numeric(table(cellchat2@idents)) # number of cells in each cell group

# Set database
cellchat2@DB <- CellChatDB

# Preprocessing
cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
# cellchat2 <- projectData(cellchat2, PPI.human)
cellchat2 <- smoothData(cellchat2, adj=PPI.human)

# Compute communication probabilities
cellchat2 <- computeCommunProb(cellchat2, raw.use = FALSE, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
df.net2 <- subsetCommunication(cellchat2)

# Compute communication probabilities at pathway level
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

# Save results
setwd(out_dir)
saveRDS(cellchat2, 'LOAD_Normal.cellchat.rds')

###################################################################
# Comparative analysis between LOAD_DEPR and LOAD_Normal
###################################################################

object.list <- list(LOAD_DEPR=cellchat1, LOAD_Normal=cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# Visualize differential number of interactions
png(filename = 'LOAD_DEPR_vs_Normal.cc.diff.number.inter.network.png', width = 2000, height = 2000, res = 300)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("LOAD_DEPR","LOAD_Normal"))
dev.off()

# Networks: red (or blue) colored edges represent increased (or decreased) signaling in LOAD_DEPR compared to LOAD_Normal

# Visualize differential interaction strength
png(filename = 'LOAD_DEPR_vs_Normal.cc.diff.inter.strength.network.png', width = 2000, height = 2000, res = 300)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c("LOAD_DEPR","LOAD_Normal"))
dev.off()

# Heatmap visualization
# In the colorbar, red (or blue) represents increased (or decreased) signaling in LOAD_DEPR compared to LOAD_Normal
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1,2))

png(filename = 'LOAD_DEPR_vs_Normal.cc.diff.number.inter.heatmap.png', width = 2000, height = 2000, res = 300)
print(gg1)
dev.off()

png(filename = 'LOAD_DEPR_vs_Normal.cc.diff.inter.strength.heatmap.png', width = 2000, height = 2000, res = 300)
print(gg2)
dev.off()

# Save merged CellChat object
saveRDS(cellchat, 'LOAD_DEPR_vs_Normal.cellchat.rds')

# Compare overall interaction numbers and strengths
gg1 <- compareInteractions(cellchat, show.legend = F, group = c("LOAD_DEPR","LOAD_Normal"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c("LOAD_DEPR","LOAD_Normal"), measure = "weight")

png(filename = 'LOAD_DEPR_vs_Normal.cc.inter.number.bar.png', width = 1000, height = 2000, res = 300)
print(gg1)
dev.off()

png(filename = 'LOAD_DEPR_vs_Normal.cc.inter.weight.bar.png', width = 1000, height = 2000, res = 300)
print(gg2)
dev.off()

# Identify signaling pathways with significant changes
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2))

png(filename = 'LOAD_DEPR_vs_Normal.cc.signal.pathways.relative.info.flow.heatmap.png', width = 1000, height = 5000, res = 200)
print(gg1)
dev.off()

png(filename = 'LOAD_DEPR_vs_Normal.cc.signal.pathways.info.flow.heatmap.png', width = 1000, height = 5000, res = 200)
print(gg2)
dev.off()

# Signaling role analysis - scatter plot showing how cell types change their roles
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

png(filename = 'LOAD_DEPR.signaling.role.scatter.png', width = 2000, height = 2000, res = 200)
print(gg[[1]])
dev.off()

png(filename = 'LOAD_Normal.signaling.role.scatter.png', width = 2000, height = 2000, res = 200)
print(gg[[2]])
dev.off()

# Create summary report
summary_file <- file.path(out_dir, "cellchat_analysis_summary.txt")
writeLines(c(
  paste("CellChat Analysis Summary - AD Depression Study"),
  paste("Date:", Sys.Date()),
  paste("Comparison: LOAD_DEPR vs LOAD_Normal"),
  paste("Dataset: 9.integrated_with_mayo_annotations.rds"),
  paste("Assay used: SCT (normalized data)"),
  paste(""),
  paste("Cell counts:"),
  paste("- Total cells after DoubletFinder filtering:", nrow(metadata)),
  paste("- Cells with clinical data:", nrow(metadata_with_clinic)),
  paste("- LOAD_DEPR cells for analysis:", length(cell.use1)),
  paste("- LOAD_Normal cells for analysis:", length(cell.use2)),
  paste(""),
  paste("Cell types with sufficient cells (>=10 in both conditions):", length(cell.type)),
  paste("Cell types analyzed:"),
  paste(paste("  -", cell.type), collapse = "\n"),
  paste(""),
  paste("Cell type counts by condition:"),
  paste(capture.output(print(cell_type_counts)), collapse = "\n"),
  paste(""),
  paste("Output files generated:"),
  paste("- LOAD_DEPR.cellchat.rds"),
  paste("- LOAD_Normal.cellchat.rds"), 
  paste("- LOAD_DEPR_vs_Normal.cellchat.rds"),
  paste("- Various visualization PNG files")
), summary_file)

print("CellChat analysis completed successfully!")
print(paste("Results saved to:", out_dir))

print(sessionInfo())