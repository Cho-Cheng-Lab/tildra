Sys.setenv(RETICULATE_PYTHON='/home/dwu/.virtualenvs/seurat/bin/python')
library(reticulate)
use_virtualenv('seurat')
py_config()
#py_install('umap-learn')
py_module_available('umap')


####

library(here)

data_dir <- here('data/derived/tildra') # data file output directory


library(Seurat) 
library(tidyverse)
library(SeuratWrappers)
library(wutilities) # devtools::install_github('symbiologist/wutilities')
library(tictoc)
library(harmony)
theme_set(theme_dwu()) # set default theme

seuratobj <- read_rds(file.path(data_dir, 'seuratobj_qc.rds'))
seuratobj

sampling <- colnames(seuratobj) %>% sample(size = 1000)
test <- subset(seuratobj, cells = sampling)

test <- test %>% 
  NormalizeData() %>%
  ScaleData() %>% 
  FindVariableFeatures(nfeatures = 500) %>% 
  RunPCA() %>% 
  FindNeighbors(dims = 1:10)

test %>% RunUMAP(dims = 1:2, umap.method  = 'umap-learn', metric = 'correlation')
test %>% RunUMAP(graph = 'RNA_snn')
test <- test %>% RunUMAP(umap.method = 'umap-learn', dims = 1:2)

#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'), plot_convergence = TRUE, assay = 'RNA')
#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'), assay = 'RNA')
#test <- test %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'))

test <- test %>% RunUMAP(reduction = 'harmony', dims = dims, ) 

p <- seurat_feature(test, features = 'RNA_snn_res.0.1')

p