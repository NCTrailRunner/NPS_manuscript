#srun --mem=100GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

set.seed(1)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(dplyr)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(glmGamPoi)
library(future)
library(harmony)
library(patchwork)

integration_method <- 'harmony'
ref <- FALSE
predicted_ids <- "4.merge_annotations___rpca_with_ref___transfer_anchors=pcaproject___wt_red=pcaproject___recomp_resids=FALSE.rds"

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/cell.type.annot.samples')
predicted_ids <- readRDS(predicted_ids)

mtx <- predicted_ids[, c('exc', 'inh', 'astro', 'micro', 'oligo', 'opc')]
x1 <- predicted_ids$prediction_score_max
x2 <- apply(mtx, MARGIN = 1, FUN = function(x){sort(x, decreasing = TRUE)[2]})
predicted_ids$grubman_score <- (x1 - x2) / x1
ix <- predicted_ids$prediction_score_max >= 0.9 & predicted_ids$grubman_score >= 0.95
predicted_ids <- predicted_ids[ix, ]

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects')
seu <- readRDS('5.sct_normalization___seurat_object.rds')
cells_use <- intersect(rownames(predicted_ids), colnames(seu))
seu <- subset(seu, cells = cells_use)
gc()

integration_fxn <- switch(integration_method, 
                          'cca' = CCAIntegration, 
                          'rpca' = RPCAIntegration, 
                          'harmony' = HarmonyIntegration)

if (ref){
  ref_ix <- sapply(seu[['RNA_raw']]@layers, 
                   function(x){
                     dim(x)[2]
                   }) %>% which.max()
} else {
  ref_ix <- NULL
}

# Create a modified batch variable
seu$batch_mod <- seu$batch

# Identify rare batches (less than 50 cells)
batch_counts <- table(seu$batch)
rare_batches <- names(batch_counts)[batch_counts < 50]
# Print the rare batches so you know what's being modified
print(rare_batches)

# Group these rare batches together
seu$batch_mod[seu$batch %in% rare_batches] <- "rare_batch"

# Run PCA first if not already done
seu <- RunPCA(seu, assay = 'SCT')

# Run Harmony directly with the modified batch variable
seu <- RunHarmony(seu, group.by.vars = "batch_mod")

#Find overlap between Seurat object and predicted_ids
shared_cells <- intersect(colnames(seu), rownames(predicted_ids))
print(paste("Cells in common:", length(shared_cells)))
print(paste("Total cells in Seurat object:", ncol(seu)))
print(paste("Total cells in predicted_ids:", nrow(predicted_ids)))

# Create subset with only matching cells
predicted_ids_subset <- predicted_ids[shared_cells, ]

# Add the predicted IDs and other metadata to the Seurat object
columns_to_add <- c("predicted_id", "prediction_score_max", "grubman_score")
seu <- AddMetaData(seu, metadata = predicted_ids_subset[, columns_to_add])

# Run UMAP on the harmony embeddings
seu <- RunUMAP(seu, dims = 1:30, reduction = 'harmony')

# Create visualization plots
# For batch plot, use the existing sample_id from seu's metadata
batch_umap <- DimPlot(seu, reduction = 'umap', group.by = 'sample_id', 
                      label = FALSE) + NoLegend()
cell_umap <- DimPlot(seu, reduction = 'umap', group.by = 'predicted_id', 
                     label = TRUE, repel = TRUE) + NoLegend()
joint_umap <- batch_umap + cell_umap

# Save the plot
setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/plots') 
png(filename = '6.inspect_umap___rna_umap.png', 
    width = 2000, height = 1000, res = 100)
print(joint_umap)
dev.off()

saveRDS(seu, paste0('6.integrate_datasets___method=', integration_method, '___ref=', ref, '.umap.rds'))

sesh <- capture.output(sessionInfo())
print(sesh)