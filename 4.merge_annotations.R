#srun --mem=100GB --pty bash -i
#conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc

.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/r_packages'))
.libPaths(c(.libPaths(), '/hpc/group/adrc/zm77/software/conda_envs'))

set.seed(1)
library(Seurat, lib.loc = '/hpc/group/adrc/zm77/r_packages')
library(dplyr)
library(Matrix, lib.loc = '/hpc/group/adrc/zm77/r_packages')

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/cell.type.annot.samples')

df <- list.files(pattern = '3.annotate_cells___.*rpca_with_ref___transfer_anchors=pcaproject___wt_red=pcaproject___recomp_resids=FALSE.rds')
df <- lapply(df, readRDS)
df <- bind_rows(df)

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/cell.type.annot.samples')
saveRDS(df, '4.merge_annotations___rpca_with_ref___transfer_anchors=pcaproject___wt_red=pcaproject___recomp_resids=FALSE.rds')

sesh <- capture.output(sessionInfo())
print(sesh)
