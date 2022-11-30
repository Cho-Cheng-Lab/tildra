# Sys.setenv(RETICULATE_PYTHON='/home/dwu/.virtualenvs/seurat/bin/python')
# library(reticulate)
# use_virtualenv('seurat')
# py_config()
# #py_install('umap-learn')
# py_module_available('umap')


####

library(here)

data_dir <- here('data/derived/tildra') # data file output directory

library(Seurat) 
library(tidyverse)
library(wutilities) # devtools::install_github('symbiologist/wutilities')
library(tictoc)
library(harmony)
theme_set(theme_dwu()) # set default theme

tic()
seuratobj <- read_rds(file.path(data_dir, 'seuratobj_qc.rds'))
seuratobj
toc()

sampling <- colnames(seuratobj) %>% sample(size = 10000)
test <- subset(seuratobj, cells = sampling)

tic()
test <- test %>% 
  NormalizeData() %>%
  ScaleData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10, return.neighbor = TRUE)
toc()

reduction <- test@reductions$pca@cell.embeddings[,1:10]

nn_object <- test@neighbors$RNA.nn

seurat_nn <- list('idx' = nn_object@nn.idx, 
                  'dist' = nn_object@nn.dist)

str(uwot)
uwot_test <- uwot::umap(X = reduction, ret_nn = TRUE)
uwot_nn <- uwot_test$nn$euclidean
str(uwot_nn)
str(seurat_nn)


uwot_nn_input <- uwot::umap(X = reduction, nn_method = seurat_nn)






test %>% RunUMAP(dims = 1:2, umap.method  = 'umap-learn', metric = 'correlation')
test %>% RunUMAP(graph = 'RNA_snn')
test <- test %>% RunUMAP(umap.method = 'umap-learn', dims = 1:2)

#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'), plot_convergence = TRUE, assay = 'RNA')
#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'), assay = 'RNA')
#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'))

test <- test %>% RunUMAP(reduction = 'harmony', dims = dims, ) 

p <- seurat_feature(test, features = 'RNA_snn_res.0.1')

p