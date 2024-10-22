## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(here)


## ----eval=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---- message=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
output_dir <- here('analysis/output/04_cluster') # analysis file output directory
data_dir <- here('data/derived/tildra') # data file output directory

dir.create(output_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)


## ---- message=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat) 
library(tidyverse)
library(ggExtra)
library(SeuratWrappers)
library(wutilities) # devtools::install_github('symbiologist/wutilities')
library(tictoc)
library(harmony)
theme_set(theme_dwu()) # set default theme


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- read_rds(file.path(data_dir, 'seuratobj_qc.rds'))
seuratobj


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- seuratobj %>% subset(subset = species == 'Human', features = str_subset(rownames(seuratobj), 'mm10---', negate = TRUE))
seuratobj


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
adt <- read_tsv(here('analysis/output/03_adt/adt_calls.tsv.gz'))

metadata_to_add <- seuratobj@meta.data %>% 
  rownames_to_column() %>% 
  select(rowname) %>% 
  left_join(adt) %>% 
  column_to_rownames('rowname')

colnames(metadata_to_add) <- paste0('ADT.', colnames(metadata_to_add))

seuratobj <- Seurat::AddMetaData(seuratobj, metadata = metadata_to_add)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- seuratobj %>% 
  NormalizeData() %>%
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- seuratobj %>% RunHarmony(group.by.vars = c('treatment', 'chemistry', 'patient', 'flow', 'mouse'), plot_convergence = TRUE, assay = 'RNA')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(here(output_dir, "elbow_plot.pdf"))
ElbowPlot(seuratobj, ndims = 50, reduction = 'harmony')
dev.off()

ElbowPlot(seuratobj, ndims = 50, reduction = 'harmony')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
map(1:20, function(i) {PCHeatmap(seuratobj, reduction = 'harmony', dims = i)})

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dims <- 1:15
neighbors <- 30
seuratobj <- seuratobj %>% 
  FindNeighbors(reduction = 'harmony', dims = dims, k.param = neighbors) %>% 
  FindClusters(resolution = c(0.1, 1)) %>% 
  RunUMAP(reduction = 'harmony', dims = dims, n.neighbors = neighbors, graph = 'RNA_snn') 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Fix cluster factor levels
clusterings <- colnames(seuratobj@meta.data) %>% str_subset('_res')

for(i in clusterings) {
  clusters <- seuratobj@meta.data[[i]]
  seuratobj[[i]] <- factor(clusters, levels = levels(clusters) %>% as.numeric() %>% sort())
}

default_clustering <- 'RNA_snn_res.0.1' 
  
seuratobj$seurat_clusters <- seuratobj[[default_clustering]]
Idents(seuratobj) <- default_clustering


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(seuratobj, features = 'RNA_snn_res.0.1', facet_hide = TRUE, legend_position = 'none', title = 'Resolution: 0.1', rasterize_dpi = 600)
ggsave(plot = p,
       filename = 'umap_res_0.1.png',
       h = 5,
       w = 5,
       dpi = 600,
       path = output_dir)
p

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(seuratobj, features = 'RNA_snn_res.0.1', facets = 'RNA_snn_res.0.1', label = FALSE, rasterize_dpi = 600)

ggsave(plot = p,
       filename = 'umap_res_0.1_facets.png',
       h = 8,
       w = 9,
       dpi = 600,
       path = output_dir)

p


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(seuratobj, 
               features = c('CD3D', #lymphocytes
                            'MS4A1', # b cells
                            'HLA-DRA', # APC
                            'LYZ', # APC
                            'TPSAB1', # mast cells
                            'COL6A2', # fibroblasts
                            'PECAM1', # endothelial cells
                            'MKI67'), # cycling
               nrow = 2, 
               rasterize_dpi = 600) 

ggsave(plot = p,
       filename = 'broad_markers.png',
       h = 4,
       w = 8,
       dpi = 600,
       path = output_dir)

p


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(seuratobj, features = 'RNA_snn_res.1', facet_hide = TRUE, color_package = 'ggplot', legend_position = 'none', title = 'Resolution: 1.0', rasterize_dpi = 600)

ggsave(plot = p,
       filename = 'umap_res_1.0.png',
       h = 5,
       w = 5,
       dpi = 600,
       path = output_dir)
p

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj %>% write_rds(here(data_dir, 'seuratobj_clustered.rds'))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

