# srun --mem=200GB --pty bash -i
# conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
# R

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs', '/hpc/group/adrc/zm77/r_packages'))

library(CellChat)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(dplyr)
options(stringsAsFactors = FALSE)	
library(NMF)
library(ggalluvial)

# Define directories
base_dir <- '/hpc/group/adrc/zm77/LUTZ.MDD/MDD_new'
out_dir <- file.path(base_dir, 'cellchat')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the integrated object with Mayo annotations
setwd(base_dir)
combined.integrated <- readRDS('objects/10.integrated_with_mayo_annotations.rds')

# Set default assay to SCT (which contains normalized data)
DefaultAssay(combined.integrated) <- "SCT"

print("Original dataset dimensions:")
print(dim(combined.integrated))

# Get metadata from Seurat object
metadata <- combined.integrated@meta.data

# Filter by DoubletFinder score if available
if("DoubletFinder.score" %in% colnames(metadata)) {
  metadata <- metadata[metadata$DoubletFinder.score < 0.75, ]
  message("After DoubletFinder filtering: ", nrow(metadata), " cells")
  
  # Update Seurat object with filtered cells
  cells_to_keep <- rownames(metadata)
  combined.integrated <- combined.integrated[, cells_to_keep]
  combined.integrated@meta.data <- metadata
} else {
  message("DoubletFinder.score not found in metadata, skipping DoubletFinder filtering")
}

# Get normalized data from SCT assay for CellChat
data.input <- GetAssayData(combined.integrated, assay = "SCT", layer = "data")

# Verify data structure
print("Data input dimensions:")
print(dim(data.input))
print("Data input class:")
print(class(data.input))

# Extract metadata
meta <- combined.integrated@meta.data

# Check condition distribution
print("Condition distribution:")
print(table(meta$Condition))

# Extract cells for each condition
cell.use1 <- rownames(meta)[meta$Condition == "Case"] # MDD cases
meta1 <- meta[cell.use1, ]
print(paste("MDD Case cells:", length(cell.use1)))
print("MDD Case cell types:")
case_cell_types <- table(meta1$mayo.cell.type, useNA = "ifany")
print(case_cell_types)

cell.use2 <- rownames(meta)[meta$Condition == "Control"] # Controls
meta2 <- meta[cell.use2, ]
print(paste("Control cells:", length(cell.use2)))
print("Control cell types:")
control_cell_types <- table(meta2$mayo.cell.type, useNA = "ifany")
print(control_cell_types)

# Find common cell types between both groups (excluding NA)
cell.type1 <- unique(meta1$mayo.cell.type)
cell.type1 <- cell.type1[!is.na(cell.type1)]
cell.type2 <- unique(meta2$mayo.cell.type) 
cell.type2 <- cell.type2[!is.na(cell.type2)]

cell.type <- intersect(cell.type1, cell.type2)
print("Common cell types between Case and Control:")
print(cell.type)

# Check cell counts for each cell type in both conditions
if(length(cell.type) > 0) {
  cell_type_counts <- data.frame(
    CellType = cell.type,
    MDD_Case = sapply(cell.type, function(x) sum(meta1$mayo.cell.type == x, na.rm = TRUE)),
    Control = sapply(cell.type, function(x) sum(meta2$mayo.cell.type == x, na.rm = TRUE))
  )
  print("Cell type counts by condition:")
  print(cell_type_counts)
  
  # Filter cell types that have at least 10 cells in both conditions (CellChat minimum)
  sufficient_cell_types <- cell_type_counts$CellType[
    cell_type_counts$MDD_Case >= 10 & cell_type_counts$Control >= 10
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
# Analysis for MDD Cases
###################################################################

# Extract cells for MDD Case condition with sufficient cell types
cell.use1 <- rownames(meta)[meta$Condition == "Case" & meta$mayo.cell.type %in% cell.type]
data.input1 <- data.input[, cell.use1]
metadata1 <- data.frame(labels = meta[cell.use1, "mayo.cell.type"], row.names = cell.use1)

print(paste("MDD Case cells for CellChat:", length(cell.use1)))
print("MDD Case cell type distribution for CellChat:")
print(table(metadata1$labels))

# Create CellChat object for MDD Cases
cellchat1 <- createCellChat(object = data.input1, meta = metadata1, group.by = "labels")
cellchat1 <- addMeta(cellchat1, meta = metadata1)
cellchat1 <- setIdent(cellchat1, ident.use = "labels")
levels(cellchat1@idents) # show factor levels of the cell labels
groupSize1 <- as.numeric(table(cellchat1@idents)) # number of cells in each cell group

# Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
str(CellChatDB)
cellchat1@DB <- CellChatDB

# Preprocessing
cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)
cellchat1 <- projectData(cellchat1, PPI.human)

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
saveRDS(cellchat1, 'MDD.Case.cellchat.rds')

###################################################################
# Analysis for Controls
###################################################################

# Extract cells for Control condition with sufficient cell types
cell.use2 <- rownames(meta)[meta$Condition == "Control" & meta$mayo.cell.type %in% cell.type]
data.input2 <- data.input[, cell.use2]
metadata2 <- data.frame(labels = meta[cell.use2, "mayo.cell.type"], row.names = cell.use2)

print(paste("Control cells for CellChat:", length(cell.use2)))
print("Control cell type distribution for CellChat:")
print(table(metadata2$labels))

# Create CellChat object for Controls
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
cellchat2 <- projectData(cellchat2, PPI.human)

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
saveRDS(cellchat2, 'MDD.Control.cellchat.rds')

###################################################################
# Comparative analysis between MDD Cases and Controls
###################################################################

object.list <- list(MDD=cellchat1, Normal=cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

# Visualize differential number of interactions
png(filename = 'MDD.cc.diff.number.inter.network.png', width = 2000, height = 2000, res = 300)
netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c("MDD","Normal"))
dev.off()

# Networks: red (or blue) colored edges represent increased (or decreased) signaling in MDD compared to Normal

# Visualize differential interaction strength
png(filename = 'MDD.cc.diff.inter.strength.network.png', width = 2000, height = 2000, res = 300)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c("MDD","Normal"))
dev.off()

# Heatmap visualization
# In the colorbar, red (or blue) represents increased (or decreased) signaling in MDD compared to Normal
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
gg2 <- netVisual_heatmap(cellchat, measure = "weight", comparison = c(1,2))

png(filename = 'MDD.cc.diff.number.inter.heatmap.png', width = 2000, height = 2000, res = 300)
print(gg1)
dev.off()

png(filename = 'MDD.cc.diff.inter.strength.heatmap.png', width = 2000, height = 2000, res = 300)
print(gg2)
dev.off()

# Save merged CellChat object
saveRDS(cellchat, 'MDD.cellchat.rds')

# Compare overall interaction numbers and strengths
gg1 <- compareInteractions(cellchat, show.legend = F, group = c("MDD","Normal"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c("MDD","Normal"), measure = "weight")

png(filename = 'MDD.cc.inter.number.bar.png', width = 1000, height = 2000, res = 300)
print(gg1)
dev.off()

png(filename = 'MDD.cc.inter.weight.bar.png', width = 1000, height = 2000, res = 300)
print(gg2)
dev.off()

# Identify signaling pathways with significant changes
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1,2))

png(filename = 'MDD.cc.signal.pathways.relative.info.flow.heatmap.png', width = 1000, height = 5000, res = 200)
print(gg1)
dev.off()

png(filename = 'MDD.cc.signal.pathways.info.flow.heatmap.png', width = 1000, height = 5000, res = 200)
print(gg2)
dev.off()

# Signaling role analysis - scatter plot showing how cell types change their roles
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets

gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

png(filename = 'MDD.signaling.role.scatter.png', width = 2000, height = 2000, res = 200)
print(gg[[1]])
dev.off()

png(filename = 'Normal.signaling.role.scatter.png', width = 2000, height = 2000, res = 200)
print(gg[[2]])
dev.off()

# Create summary report
summary_file <- file.path(out_dir, "cellchat_analysis_summary.txt")
writeLines(c(
  paste("CellChat Analysis Summary - MDD Study"),
  paste("Date:", Sys.Date()),
  paste("Comparison: MDD Cases vs Controls"),
  paste("Dataset: 10.integrated_with_mayo_annotations.rds"),
  paste("Assay used: SCT (normalized data)"),
  paste(""),
  paste("Cell counts:"),
  paste("- Total cells in dataset:", ncol(combined.integrated)),
  paste("- MDD Case cells for analysis:", length(cell.use1)),
  paste("- Control cells for analysis:", length(cell.use2)),
  paste(""),
  paste("Cell types with sufficient cells (>=10 in both conditions):", length(cell.type)),
  paste("Cell types analyzed:"),
  paste(paste("  -", cell.type), collapse = "\n"),
  paste(""),
  paste("Cell type counts by condition:"),
  paste(capture.output(print(cell_type_counts)), collapse = "\n"),
  paste(""),
  paste("Output files generated:"),
  paste("- MDD.Case.cellchat.rds"),
  paste("- MDD.Control.cellchat.rds"), 
  paste("- MDD.cellchat.rds"),
  paste("- Various visualization PNG files")
), summary_file)

print("CellChat analysis completed successfully!")
print(paste("Results saved to:", out_dir))

print(sessionInfo())