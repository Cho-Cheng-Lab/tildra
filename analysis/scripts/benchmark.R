library(tictoc)

tic()
library(here)

data_dir <- here('data/derived/tildra') # data file output directory

library(Seurat) 
library(tidyverse)
library(wutilities) # devtools::install_github('symbiologist/wutilities')
library(tictoc)
library(harmony)
theme_set(theme_dwu()) # set default theme

tic()
seuratobj <- read_rds(file.path(data_dir, 'seuratobj_clustered.rds'))
seuratobj
toc()

sampling <- colnames(seuratobj) %>% sample(size = 150000)
test <- subset(seuratobj, cells = sampling)

tic()
test <- test %>% 
  NormalizeData() %>%
  ScaleData() %>% 
  FindVariableFeatures(nfeatures = 1000) %>% 
  RunPCA() 

toc()

tic()
test <- test %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters() %>% 
  RunUMAP(reduction = 'pca',
          dims = 1:10)
toc()

toc()