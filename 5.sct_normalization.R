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

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects/scSampler/downsampled_samples')
seu <- list.files(pattern = '.rds')

seu <- lapply(seu, function(x){
  obj <- readRDS(x)
  DefaultAssay(obj) <- 'RNA_raw'
  return(obj)
})

ix <- sapply(seu, ncol) >= 50 # sample R1583702 not enough cells
seu <- seu[ix]
seu <- merge(seu[[1]], seu[2:length(seu)])

options(future.globals.maxSize = 20 * 1024^3) 
seu <- SCTransform(seu, vst.flavor = 'v2', 
                   variable.features.n = 5000, assay = 'RNA_raw')

setwd('/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/objects')
saveRDS(seu, '5.sct_normalization___seurat_object.rds')


sesh <- capture.output(sessionInfo())
print(sesh)