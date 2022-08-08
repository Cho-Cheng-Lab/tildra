library(Seurat) 
library(tidyverse)
library(wutilities)
library(tictoc)
library(parallel)

output_dir <- 'analysis/output/03_markers' # analysis file output directory
data_dir <- 'data/derived/seurat' # data file output directory

seuratobj <- read_rds(file.path(data_dir, 'seuratobj_clustered.rds'))
markers <- read_tsv(file.path(output_dir, 'markers', 'markers_res_0.2.tsv'))

genes_to_plot <- markers %>% filter(avg_log2FC > 1, !str_detect(gene, 'mm10---')) %>% pull(gene) %>% unique()
completed_genes <- list.files(plot_dir) %>% str_remove('.png')
remaining_genes <- setdiff(genes_to_plot, completed_genes)

tic()
mclapply(1:length(remaining_genes),
         mc.cores = 40,
         FUN = function(i) {
           
           gene <- remaining_genes[i]
           print2(paste0('Plotting gene ', i, ' of ', length(remaining_genes), ': ', gene))
           p <- seurat_feature(seuratobj, 
                               features = gene)
           
           ggsave(plot = p,
                  filename = paste0(gene, '.png'),
                  path = plot_dir,
                  h = 3,
                  w = 3,
                  dpi = 75)
           
         })
toc()