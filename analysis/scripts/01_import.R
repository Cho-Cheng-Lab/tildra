## -------------------------------------------------------------------------------------------------------------------------
library(here)


## ----eval=FALSE, message=FALSE--------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---- message=FALSE-------------------------------------------------------------------------------------------------------
output_dir <- here('analysis/output/01_import') # analysis file output directory
data_dir <- here('data/derived/tildra') # data file output directory

dir.create(output_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)


## ---- message=FALSE-------------------------------------------------------------------------------------------------------
library(Seurat) 
library(tidyverse)
library(scDblFinder) 
library(patchwork)
library(ggExtra)
library(ggplotify)
library(parallel)
library(tictoc)
library(wutilities) # devtools::install_github('symbiologist/wutilities')


## -------------------------------------------------------------------------------------------------------------------------
#function to trim ADT:
trim.feature.names <- function(inmat){
  
  rownames(inmat) <- rownames(inmat) %>% str_remove('_TotalA') %>% str_remove('/')
  
  return(inmat)
}

#qc plot
qc_plots <- function(seuratobj,
                     filename,
                     path) {
  
  dir.create(path = path, showWarnings = FALSE, recursive = TRUE)
  
  p1 <-  seuratobj@meta.data %>% 
    filter(pct_mouse < 5) %>% 
    ggplot(aes(x = nCount_RNA, 
               y = nFeature_RNA)) +
    ggtitle(filename) + 
    geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
    geom_hline(yintercept=2500, 
               linetype="dashed", 
               color = "black", 
               size=1) + 
    geom_vline(xintercept = 10000, 
               linetype = "dashed", 
               color = "black", 
               size=1) +
    theme_bw()
  
  p1 <- ggMarginal(p1, type = 'densigram')
  
  p2 <- seuratobj@meta.data %>% 
    filter(pct_mouse < 5) %>% 
    ggplot(aes(x = pct_ribo, 
               y = pct_mito)) +
    ggtitle(filename) + 
    geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
    geom_hline(yintercept = 20, linetype = "dashed", color = "black", size = 1) + 
    geom_vline(xintercept = 50, linetype="dashed", color = "black", size = 1) +
    theme_bw()
    
  p2 <- ggMarginal(p2, type = 'densigram')

    p3 <- seuratobj@meta.data %>% 
    filter(pct_mouse < 5) %>% 
    ggplot(aes(x = pct_stress, 
               y = pct_mito)) +
    ggtitle(filename) + 
    geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
    geom_hline(yintercept = 20, linetype = "dashed", color = "black", size = 1) + 
    geom_vline(xintercept = 30, linetype="dashed", color = "black", size = 1) +
    theme_bw()
    
        species_tally <- seuratobj@meta.data %>% 
      mutate(species = factor(ifelse(pct_human > 95, 'Human', 'Mouse'), levels = c('Human', 'Mouse'))) %>% 
      group_by(species, .drop = FALSE) %>% 
      tally() %>% 
      mutate(pct = round(100 * n / sum(n), 2))
    
    doublet_tally <- seuratobj@meta.data %>% 
      group_by(doublet_class) %>% 
      tally() %>% 
      mutate(pct = round(100 * n / sum(n), 2))
    
    p4 <- species_tally %>% 
      ggplot(aes(x = species,
                 y = pct,
                 label = paste0('n = ',n, ', ', pct, '%'),
                 fill = species)) +
      geom_col(color = 'black') +
      geom_label(fill = 'white') +
      ggthemes::scale_fill_few() +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(title = 'Species Assignment', 
           y = 'Percent of Cells',
           x = 'Species') +
      ylim(c(0, 100))
    
    p5 <- doublet_tally %>% 
      ggplot(aes(x = doublet_class,
                 y = pct,
                 label = paste0('n = ',n, ', ', pct, '%'),
                 fill = doublet_class)) +
      geom_col(color = 'black') +
      geom_label(fill = 'white') +
      ggthemes::scale_fill_few(palette = 'Light') +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(title = 'Doublet Assignment',
           y = 'Percent of Cells',
           x = 'Doublet Status') +
      ylim(c(0, 100))
      
    ### histograph
    page1 <- list(p1, p2) %>% wrap_plots(nrow = 1)
    page2 <- p3 | (p4 + p5)
    
    pdf(paste0(path, '/', filename, '.pdf'), width=14)
    options(repr.plot.height = 6, repr.plot.width = 15)
  
    print(page1)
    print(page2)
    
    options(repr.plot.height = 6, repr.plot.width = 24)
    plot1 <-  hist(seuratobj$nCount_RNA, breaks = 100, main=filename, xlab = "Count depth (Counts/cell)", col='wheat')
    plot1 <-  abline(v=400,col="red")
    plot2 <-  hist(seuratobj$nFeature_RNA, breaks = 200, main=filename, xlab = "Number of genes per cell", xlim = c(0, 4000), col='wheat')
    plot2 <-  abline(v=250,col="red")
    print(plot1 + plot2)
    
    dev.off()
  
}


## -------------------------------------------------------------------------------------------------------------------------
samples_dir <- here('analysis/input/samples/')
metadata <- read_tsv(here('analysis/output/00_process/metadata_subset.tsv'))
stress_genes <- read_csv(here('analysis/input/coregene_df-FALSE-v3.csv'))


## -------------------------------------------------------------------------------------------------------------------------
samples <- metadata$id %>% set_names(.) 

tic()
seurats <- 
  mclapply(samples, mc.cores = 10, # in parallel
  #map(samples,         
      function(i) {
        
        print2(i)
        
        # Sample metadata
        sample_subset <- metadata %>% filter(id == i) %>% as.data.frame()
        metadata_to_add <- sample_subset %>% 
          select(id,
                 sample,
                 patient,
                 condition,
                 treatment,
                 response,
                 group,
                 sample_label,
                 patient_label, 
                 adt,
                 flow,
                 mouse,
                 saturation,
                 chemistry,
                 contains('ab_'),
                 published,
                 sample_order)
        
        # read in data from 10x output
        tenx_matrix <- Read10X_h5(here(sample_subset$h5))
        
        if(length(tenx_matrix) == 2) { # two assays
          # create seurat obj
          seuratobj <- CreateSeuratObject(counts = CollapseSpeciesExpressionMatrix(tenx_matrix[['Gene Expression']], prefix = 'GRCh38_', controls = 'mm10___', ncontrols = 200), project = i)
          
          # add additional ADT assay
          seuratobj[['ADT']] <- CreateAssayObject(counts = trim.feature.names(tenx_matrix[['Antibody Capture']]))
          
        } else { # single assay, gene expression only
          seuratobj <- CreateSeuratObject(counts = tenx_matrix,
                                          project = i)
        }
        
        # filter out low complexity
        seuratobj <- subset(seuratobj, subset = nFeature_RNA > 100)
          
        # calculate doublet class
        sceobj <- as.SingleCellExperiment(seuratobj, assay = 'RNA')
        sceobj <- scDblFinder(sceobj, verbose = FALSE)
        
        # add metadata
        seuratobj <- AddMetaData(seuratobj, metadata = rep(metadata_to_add, ncol(seuratobj)))
        
        # add QC metrics
        seuratobj <- PercentageFeatureSet(seuratobj, pattern = '^MT-|mt-', col.name = 'pct_mito')
        seuratobj <- PercentageFeatureSet(seuratobj, pattern = '^RP[SL]|^MRP[SL]|Rp[sl]|Mrp[sl]', col.name = 'pct_ribo')
        seuratobj <- PercentageFeatureSet(seuratobj, features = intersect(stress_genes$gene_symbol, rownames(seuratobj)), col.name = 'pct_stress')
        seuratobj <- PercentageFeatureSet(seuratobj, pattern = '^mm10', col.name = 'pct_mouse')
        seuratobj$pct_human <- round(100 - seuratobj$pct_mouse, 2)
        seuratobj$species <- case_when(
          seuratobj$pct_human > 95 ~ 'Human',
          seuratobj$pct_mouse > 95 ~ 'Mouse',
          TRUE ~ 'Mixed'
        )
        seuratobj$doublet_class <- sceobj$scDblFinder.class
        
        # plot QC metrics
        qc_plots(seuratobj = seuratobj,
                 filename = i,
                 path = here(output_dir, 'qc_plots'))
        
        # return seurat object
        seuratobj
        
      })
toc()


## -------------------------------------------------------------------------------------------------------------------------
seurats_merged <- merge(seurats[[1]], seurats[-1])
seurats_merged

## -------------------------------------------------------------------------------------------------------------------------
qc_table <- seurats_merged@meta.data %>% 
  filter(species == 'Human') %>% 
  select('id', ends_with('RNA'), starts_with('pct')) %>% 
  pivot_longer(cols = -1,
               names_to = 'metric',
               values_to = 'value')
  
qc_table %>% head()


## ----fig.width = 12, fig.height=16----------------------------------------------------------------------------------------
qc_plots_pre <- map(qc_table$metric %>% unique(), function(i){
  
  print2(i)
  p <- qc_table %>% 
    filter(metric == i) %>% 
    ggplot(aes(x = id,
               y = value,
               fill = id)) +
    geom_violin() +
    theme_dwu() +
    labs(x = '',
         y = '',
         title = i)
  
  if(str_detect(i, 'RNA')) {
   p <- p + scale_y_log10()
  }
  
  p
  
}) %>% wrap_plots(ncol = 1)

ggsave(plot = qc_plots_pre,
       path = output_dir,
       filename = 'qc_summary_pre.png',
       h = 16,
       w = 12,
       dpi = 150)

qc_plots_pre


## -------------------------------------------------------------------------------------------------------------------------
seurats_merged %>% write_rds(here(data_dir, 'seurats_merged.rds'))


## -------------------------------------------------------------------------------------------------------------------------
sessionInfo()

