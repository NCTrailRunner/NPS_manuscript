#srun --mem=200GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

set.seed(1)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(dplyr)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')

i <- Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects')
ref <- "4.inspect_umap___rpca_with_ref_reference_data.rds"
recomp_resids <- FALSE
ann_method <- 'pcaproject'
wt_red <- 'pcaproject'
int_method <- gsub('.*___|_reference_data.rds', '', ref)
ref <- readRDS(ref)

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/scSampler/downsampled_samples')
file_list <- list.files(pattern = '.rds')

que <- file_list[i]
sample_id <- gsub('_downsampled.rds', '', que)
message(i, ': ', sample_id)

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/scSampler/downsampled_samples')
que <- readRDS(que)
#que <- RenameAssays(que, assay.name = 'RNA_raw', new.assay.name = 'RNA')
DefaultAssay(que) <- 'RNA_raw'
que <- SCTransform(que, vst.flavor = 'v2', 
				 variable.features.n = 5000, assay = 'RNA_raw')
anchor_features <- rownames(ref[['integrated']]@scale.data)
mtx <- FetchResiduals(que, features = anchor_features, umi.assay = 'RNA_raw')
anchor_features <- rownames(que[['SCT']]@data)
anchor_features <- anchor_features[anchor_features %in% rownames(mtx)]
que[['SCT']]@scale.data <- mtx[anchor_features, ]

que <- RunPCA(que, features = anchor_features, assay = 'SCT')

que$log_umi <- log10(que$nCount_RNA_raw) # if leave out, causes error for some samples.

transfer_anchors <- FindTransferAnchors(reference = ref, query = que, 
									  normalization.method = 'SCT', 
									  recompute.residuals = recomp_resids, 
									  reference.assay = 'integrated', features = anchor_features,
									  query.assay = 'SCT', reduction = ann_method)

n <- length(unique(transfer_anchors@anchors[, 'cell2']))
k <- min( c(50, n-1) )

que <- TransferData(anchorset = transfer_anchors, query = que,
				  refdata = ref$cell_type, weight.reduction = wt_red,
				  k.weight = k)
predicted_id <- que[['prediction.score.id']]@data
predicted_id <- t(predicted_id)
predicted_id <- as.data.frame(predicted_id)
ix <- apply(predicted_id, MARGIN = 1, FUN = which.max)
predicted_id$predicted_id <- colnames(predicted_id)[ix]
ix <- cbind(1:nrow(predicted_id), ix)
predicted_id$prediction_score_max <- as.numeric(predicted_id[ix])

file_name <- paste0('3.annotate_cells___', sample_id, 
				  '___ref=', int_method, 
				  '___transfer_anchors=', ann_method, 
				  '___wt_red=', wt_red, 
				  '___recomp_resids=', recomp_resids, '.rds')
setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/cell.type.annot.samples')
saveRDS(object = predicted_id, file = file_name)

sesh <- capture.output(sessionInfo())
print(sesh)