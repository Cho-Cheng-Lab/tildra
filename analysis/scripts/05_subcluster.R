## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(here)


## ----eval=FALSE, message=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---- message=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
output_dir <- here('analysis/output/05_subcluster') # analysis file output directory
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


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
recluster <- function(seuratobj, 
                      assay = 'RNA',
                      dims = 1:20,
                      resolution = c(0.2, 1.0, 2.0),
                      default_ident = 'RNA_snn_res.1') {
  
    ### RunUMAP, FindNeighbors, and FindClusters 
  
  seuratobj <- seuratobj %>% 
    RunUMAP(reduction = "harmony", dims = dims, assay = assay) %>% 
    FindNeighbors(reduction = "harmony", dims = dims, assay = assay) %>% 
    FindClusters(resolution = resolution)
  
  ### Fix cluster factor levels
  clusterings <- colnames(seuratobj@meta.data) %>% str_subset('_res')
  
  for(i in clusterings) {
    clusters <- seuratobj@meta.data[[i]]
    seuratobj[[i]] <- factor(clusters, levels = levels(clusters) %>% as.numeric() %>% sort())
  }
  
  Idents(seuratobj) <- default_ident
  seuratobj$seurat_clusters <- seuratobj[[default_ident]]
  
  seuratobj
}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- read_rds(file.path(data_dir, 'seuratobj_clustered.rds'))
seuratobj


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(seuratobj, features = c('RNA_snn_res.0.2'), facet_hide = TRUE, legend_position = 'none', title = 'Resolution: 0.2')

p


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_feature(seuratobj, feature = 'RNA_snn_res.0.2', facets = 'RNA_snn_res.0.2', label = FALSE)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
markers_reynolds <- c('KRT14', 'KRT5', 
                      'GATA3', 'KRTDAP',
                      'CRYAB', 'TYRP1',
                      'DKK3', 'PLEKHA4',
                      'LUM', 'APOD',
                      'PHLDA2', 'CCL2',
                      'IL6', 'PDGFRA',
                      'COL1A1', 'CXCL12',
                      'CD82', 'CCL19', 'KCNE4', 'CPE', 'TAGLN', 'MYL9', 'HEY1', 'CCL14', 'PECAM1',
                      'ACKR1', 'SELE', 'SNCG', 'HES1', 'MMRN1', 'XCL1' ,'XCL2', 'KLRC1', 'CCL5', 'GZMA', 
                      'CD8A', 'CD8B', 'IL7R', 'FOXP3', 'TIGIT', 'TPSB2', 'TPSAB1', 'CD79A', 'JCHAIN', 'CD163',
                      'C1QB', 'FCGR2A', 'MS4A6A', 'IL23A', 'CLEC9A', 'CLEC10A', 'CD1C', 'CD207', 'CD1A', 
                      'CD14', 'IL1B', 'CCR7', 'IDO1')




## ----fig.height=10, fig.width=10---------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(seuratobj) <- 'RNA_snn_res.0.2'
p <- DotPlot(seuratobj, features = rev(markers_reynolds), cluster.idents = FALSE) + coord_flip()
p


## ----fig.height = 3, fig.width = 12------------------------------------------------------------------------------------------------------------------------------------------------------------
markers_keratinocytes <- c('KRT5', 'KRT14', 'GATA3', 'KRTDAP')
seurat_feature(seuratobj, features = markers_keratinocytes, nrow = 1)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_feature(seuratobj, 
               features = c('CD3D', #lymphocytes
                            'MS4A1', # b cells
                            'TPSAB1', # mast cells
                            'COL6A2'),
               nrow = 1) # endothelial cells

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_feature(seuratobj, features = 'RNA_snn_res.0.2')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seuratobj <- seuratobj %>% 
  FindNeighbors(reduction = 'umap', 
                dims = 1:2) %>% 
  FindClusters(resolution = 0.01)

seuratobj %>% seurat_feature(features = 'seurat_clusters', facet_hide = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
broad_markers <- RunPrestoAll(seuratobj) %>% mutate(pct.diff = pct.1 - pct.2) %>% arrange(-pct.diff) %>% as_tibble() %>% select(cluster, gene, pct.1, pct.2, pct.diff, everything())
broad_markers

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
broad_markers %>% group_by(cluster) %>% top_n(5, pct.diff) %>% arrange(cluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
batch <- c('treatment', 'patient', 'flow', 'mouse', 'chemistry') # global batch variables

apc_clusters <- 1
lymphocyte_clusters <- c(0, 2, 3)

apcs <- subset(seuratobj, subset = seurat_clusters %in% apc_clusters) 
lymphocytes <- subset(seuratobj, subset = seurat_clusters %in% lymphocyte_clusters)

p1 <- seurat_feature(apcs, facet_hide = TRUE, title = 'APCs')
p2 <- seurat_feature(lymphocytes, facet_hide = TRUE, title = 'Lymphocytes')
p1 + p2


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
apcs <- apcs %>% DietSeurat() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA() %>% RunHarmony(batch)
apcs


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ElbowPlot(apcs, ndims = 50, reduction = 'harmony')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lymphocytes <- lymphocytes %>% DietSeurat() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA() %>% RunHarmony(batch)
lymphocytes


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ElbowPlot(lymphocytes, ndims = 50, reduction = 'harmony')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dims <- 1:20
apcs <- apcs %>% recluster(dims = dims)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(apcs, facet_hide = TRUE, title = 'APC Subcluster', legend_position = 'none')

ggsave(plot = p,
       filename = 'umap_apc.png',
       h = 5,
       w = 5,
       path = output_dir)
p


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
apcs@meta.data %>% 
  select(seurat_clusters, id) %>% 
  unique() %>% 
  add_count(seurat_clusters) %>% 
  arrange(n)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lymphocytes <- lymphocytes %>% recluster(dims = dims)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- seurat_feature(lymphocytes, facet_hide = TRUE)
p

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lymphocytes@meta.data %>%
  select(seurat_clusters, id) %>%
  unique() %>%
  add_count(seurat_clusters) %>%
  arrange(n)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seurat_feature(lymphocytes, features = c('JCHAIN', 'MS4A1'))

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
VlnPlot(lymphocytes, 'MS4A1')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster <- subset(seuratobj, cells = c(colnames(apcs), colnames(lymphocytes))) %>% DietSeurat()

combined_features <- union(apcs@assays$RNA@var.features, lymphocytes@assays$RNA@var.features)
VariableFeatures(subcluster) <- combined_features

subcluster <- subcluster %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = batch)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
ElbowPlot(subcluster, ndims = 50, reduction = 'harmony')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster <- subcluster %>% recluster(dims = 1:20,
                                       resolution = c(0.2, 1, 2),
                                       default_ident = 'RNA_snn_res.2')

subcluster$global <- Idents(subcluster)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster %>% seurat_feature()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster@meta.data %>% 
  select(cluster = RNA_snn_res.2, id, patient, sample) %>% 
  unique() %>% 
  add_count(cluster) %>% 
  arrange(n)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
apc_clustering <- apcs@meta.data %>% 
  rownames_to_column() %>% 
  select(rowname, seurat_clusters) %>% 
  mutate(subcluster = factor(paste0('A', seurat_clusters), levels = paste0('A', levels(seurat_clusters))))

lymphocyte_clustering <- lymphocytes@meta.data %>% 
  rownames_to_column() %>% 
  select(rowname, seurat_clusters) %>% 
  mutate(subcluster = factor(paste0('L', seurat_clusters), levels = paste0('L', levels(seurat_clusters))))

merged_clustering <- subcluster@meta.data %>% 
  rownames_to_column() %>% 
  select(rowname) %>% left_join(
    bind_rows(lymphocyte_clustering,
              apc_clustering) %>% 
      select(-seurat_clusters)) %>% 
  column_to_rownames()

merged_clustering

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster <- AddMetaData(subcluster, merged_clustering)
subcluster$seurat_clusters <- subcluster$subcluster
Idents(subcluster) <- 'seurat_clusters'



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample_cluster_tally <- subcluster@meta.data %>% 
  select(cluster = RNA_snn_res.2, id, patient, sample) %>% 
  unique() %>% 
  add_count(cluster) %>% 
  arrange(n)

sample_cluster_tally


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
clusters_to_remove <- sample_cluster_tally %>% filter(n < 3) %>% pull(cluster) %>% as.character()
clusters_to_remove

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sample_subcluster_tally <- subcluster@meta.data %>% 
  filter(!(RNA_snn_res.2 %in% clusters_to_remove)) %>% 
  select(cluster = seurat_clusters, id, patient, sample) %>% 
  unique() %>% 
  add_count(cluster) %>% 
  arrange(n)

sample_subcluster_tally

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subclusters_to_remove <- sample_subcluster_tally %>% filter(n < 3) %>% pull(cluster) %>% as.character()
subclusters_to_remove


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cells_to_keep <- subcluster@meta.data %>% 
  filter(!(RNA_snn_res.2 %in% clusters_to_remove),
         !(subcluster %in% subclusters_to_remove)) %>% 
  rownames()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster <- subset(subcluster, cells = cells_to_keep)
subcluster


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster %>% seurat_feature(feature = 'subcluster', facet_hide = TRUE, title = 'Lymphocyte and APC Clusters', legend_position = 'none')


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p3 <- subcluster %>% seurat_feature(facet_hide = TRUE, title = 'Lymphocyte and APC clusters', legend_position = 'none', color_package = 'carto', color_palette = 'Bold', alpha = 0.2, label_size = 4)

ggsave(plot = p3,
       filename = 'umap_unified_labeled.png',
       h = 5,
       w = 5,
       path = output_dir)

p3

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p4 <- subcluster %>% seurat_feature(facet_hide = TRUE, title = 'Lymphocyte and APC clusters', legend_position = 'none', color_package = 'carto', color_palette = 'Bold', alpha = 0.2, label = FALSE)

ggsave(plot = p4,
       filename = 'umap_unified.png',
       h = 5,
       w = 5,
       path = output_dir)

p4


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p5 <- subcluster %>% seurat_feature(title = 'Lymphocyte and APC clusters', facets = 'seurat_clusters', 
                                    legend_position = 'none', color_package = 'carto', color_palette = 'Bold', alpha = 0.2, label = FALSE)

ggsave(plot = p5,
       filename = 'umap_unified_facet_subclusters.png',
       h = 12,
       w = 12,
       path = output_dir)

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p6 <- subcluster %>% seurat_feature(title = 'Lymphocyte and APC clusters', facets = 'RNA_snn_res.2', 
                                    legend_position = 'none', color_package = 'carto', color_palette = 'Pastel', alpha = 0.2, label = FALSE)

ggsave(plot = p6,
       filename = 'umap_unified_facet_clusters.png',
       h = 12,
       w = 12,
       path = output_dir)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster %>% write_rds(here(data_dir, 'seuratobj_subcluster.rds'))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subcluster %>% seurat_feature(facet_hide = TRUE, title = 'Lymphocyte and APC clusters', legend_position = 'none', color_package = 'carto', color_palette = 'Bold', alpha = 0.2, label_size = 4)

