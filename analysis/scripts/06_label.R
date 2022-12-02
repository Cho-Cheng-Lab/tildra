## -------------------------------------------------------------------------------------------------------------------------
library(here)


## ----eval=FALSE, message=FALSE--------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---- message=FALSE-------------------------------------------------------------------------------------------------------
output_dir <- 'analysis/output/06_label' # analysis file output directory
data_dir <- 'data/derived/tildra' # data file output directory

dir.create(output_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)


## ---- message=FALSE-------------------------------------------------------------------------------------------------------
library(Seurat) 
library(tidyverse)
library(wutilities) # devtools::install_github('symbiologist/wutilities')
library(tictoc)
theme_set(wutilities::theme_dwu())


## -------------------------------------------------------------------------------------------------------------------------
seuratobj <- read_rds(here(data_dir, 'seuratobj_subcluster.rds'))
seuratobj


## -------------------------------------------------------------------------------------------------------------------------
rashx <- read_rds(here('data/external/rashx.rds'))
rashx


## -------------------------------------------------------------------------------------------------------------------------
library(scibet)


## -------------------------------------------------------------------------------------------------------------------------
shared_genes <- intersect(rownames(seuratobj), rownames(rashx))

tic()
ref_input <- rashx@assays$RNA@data[shared_genes, ] %>% as.matrix() %>% t() %>% as.data.frame()
ref_input <- cbind(ref_input, 'label' = rashx$Ident2)
query_input <- seuratobj@assays$RNA@data[shared_genes, ] %>% as.matrix() %>% t() %>% as.data.frame()
toc()


## -------------------------------------------------------------------------------------------------------------------------
tic()
prd <- SciBet(train = ref_input,
              test = query_input, 
              k = 3000)

seuratobj$label <- prd
toc()


## -------------------------------------------------------------------------------------------------------------------------
p1 <- seurat_feature(seuratobj, features = 'label', facet_hide = TRUE, color_package = 'carto', color_palette = 'Bold', rasterize_dpi = 600, legend_position = 'none')

ggsave(plot = p1,
       filename = 'umap_label.png',
       path = output_dir,
       dpi = 600,
       w = 5, 
       h = 5)

p1


## -------------------------------------------------------------------------------------------------------------------------
p2 <- seurat_feature(seuratobj, features = 'local', facets = 'label', facet_hide = FALSE, color_package = 'carto', color_palette = 'Bold', label = FALSE, legend_position = 'none', rasterize_dpi = 600)

ggsave(plot = p2,
       filename = 'umap_label_facets.png',
       path = output_dir,
       dpi = 600,
       w = 10, 
       h = 10)
p2


## -------------------------------------------------------------------------------------------------------------------------
seuratobj@meta.data %>% rownames_to_column() %>% select(rowname, label) %>% write_tsv(here(output_dir, 'labels.tsv.gz'))


## -------------------------------------------------------------------------------------------------------------------------
rm(ref_input)
rm(query_input)
gc()


## -------------------------------------------------------------------------------------------------------------------------
library(future)
options(future.globals.maxSize = 60 * 1024^3)
plan('multicore', workers = 60)
 
transfer_dims <- 1:30
tic()
anchors <- FindTransferAnchors(reference = rashx, 
                               query = seuratobj,
                               dims = transfer_dims, 
                               features = VariableFeatures(seuratobj),
                               reduction = 'pcaproject',
                               reference.reduction = 'pca')


transfer <- TransferData(anchorset = anchors,
                            refdata = rashx$Ident2,
                            dims = transfer_dims)
toc()

transfer %>% rownames_to_column() %>% write_tsv(here(output_dir, 'transfer.tsv.gz'))


## -------------------------------------------------------------------------------------------------------------------------
transfer

## -------------------------------------------------------------------------------------------------------------------------
seuratobj$transfer <- transfer$predicted.id


## -------------------------------------------------------------------------------------------------------------------------
p1 <- seurat_feature(seuratobj, features = 'transfer', facet_hide = TRUE, color_package = 'carto', color_palette = 'Bold', rasterize_dpi = 600, legend_position = 'none')

ggsave(plot = p1,
       filename = 'umap_transfer.png',
       path = output_dir,
       dpi = 600,
       w = 5, 
       h = 5)

p1


## -------------------------------------------------------------------------------------------------------------------------
p2 <- seurat_feature(seuratobj, features = 'local', facets = 'transfer', facet_hide = FALSE, color_package = 'carto', color_palette = 'Bold', label = FALSE, legend_position = 'none', rasterize_dpi = 600)

ggsave(plot = p2,
       filename = 'umap_transfer_facets.png',
       path = output_dir,
       dpi = 600,
       w = 10, 
       h = 10)
p2


## -------------------------------------------------------------------------------------------------------------------------
sessionInfo()

