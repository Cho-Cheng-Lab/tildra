## ----setup---------------------------------------------------------------------------------------------------------------------------------------------------------------------
require(knitr)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())


## ----eval=FALSE, message=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------
## current_file <- rstudioapi::getActiveDocumentContext()$path
## output_file <- stringr::str_replace(current_file, '.Rmd', '.R')
## knitr::purl(current_file, output = output_file)
## file.edit(output_file)


## ---- message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------
output_dir <- 'analysis/output/01_import' # analysis file output directory
data_dir <- 'data/derived/seurat' # data file output directory

dir.create(output_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)


## ---- message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat) 
library(tidyverse)
library(scDblFinder) 
library(ggExtra)
library(ggplotify)


## ---- message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------
# library(metap) 
# library(MAST) 
# 
# library(ggplot2) 
# library(cowplot) 
# library(patchwork) 
# library(writexl)
# library(dplyr)
# library(xlsx)
# library(harmony)
# library(openxlsx)
# library(rgl)
# library(ggExtra)
# library(ggplotify)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
samples_dir <- "analysis/input/samples/"
metadata_df <- read_tsv('analysis/output/00_process/metadata_subset.tsv')
stress_genes <- read_csv('analysis/input/coregene_df-FALSE-v3.csv')
dir.create(file.path(output_dir, 'qc_plots', 'scatter_plot_per_sample'), recursive = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#setwd(samples_dir)
all_samples <- list()
all_samples_data <- list()
all_samples_sce <- list()

#function to trim ADT:
trim.feature.names <- function(inmat){
  newnames <- sapply(strsplit(rownames(inmat), "_"),
                     function(x) {
                       if(length(x)==1) return(x)
                       else return(paste(x[-length(x)], collapse="_"))
                     })
  rownames(inmat) <- newnames
  return(inmat)
}


## ----echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------
for(f in metadata_df$id){
	print(f)

	# read in the samples and remove doublets
	all_samples_data[[f]] <- Read10X_h5(paste0(samples_dir, f, '.h5'))

	if(!is.null(names(all_samples_data[[f]]))){
		all_samples_data[[f]][["Gene Expression"]] <- CollapseSpeciesExpressionMatrix(all_samples_data[[f]][["Gene Expression"]], prefix = "GRCh38_", controls = "mm10___", ncontrols = 500)
		all_samples_data[[f]][["Antibody Capture"]] <- trim.feature.names(all_samples_data[[f]][["Antibody Capture"]])	
		all_samples[[f]] <- CreateSeuratObject(counts = all_samples_data[[f]][["Gene Expression"]], project = as.character(f))
	} else {
		all_samples[[f]] <- CreateSeuratObject(counts = all_samples_data[[f]])
	}

	all_samples_sce[[f]] <- as.SingleCellExperiment(all_samples[[f]])
	all_samples_sce[[f]] <- scDblFinder(all_samples_sce[[f]])
	print(table(call=all_samples_sce[[f]]$scDblFinder.class))

	all_samples[[f]] <- as.Seurat(all_samples_sce[[f]], project = f)	
	all_samples[[f]][["percent.mt"]] <- PercentageFeatureSet(all_samples[[f]], pattern = "^MT-|mt-", assay = "RNA")
	all_samples[[f]][["percent.ribo"]] <- PercentageFeatureSet(all_samples[[f]], pattern = "^RP[SL]|^MRP[SL]|Rp[sl]|Mrp[sl]", assay = "RNA")
	all_samples[[f]][["percent.stress"]] <- PercentageFeatureSet(all_samples[[f]], features = intersect(stress_genes$gene_symbol, rownames(all_samples[[f]])), assay="RNA")
	all_samples[[f]][["percent.mouse"]] <- PercentageFeatureSet(all_samples[[f]], pattern = '^mm10', assay = "RNA")
	
	if(!is.null(names(all_samples_data[[f]]))){
		all_samples[[f]][["ADT"]] <- CreateAssayObject(counts = all_samples_data[[f]][["Antibody Capture"]])
	}

	all_samples[[f]] <- subset(all_samples[[f]], scDblFinder.class == "singlet")
	
	# add the metadata
	tmp_df <- subset(metadata_df, id == str_remove(f, '.h5'))
	print(paste0("adding metadata for sample ",f))

	for(m in colnames(tmp_df)){
		print(paste0("adding ",m))
		all_samples[[f]][[m]] <- tmp_df[[m]]
	}
	
	## scatterplot
	pdf(paste0(output_dir,"/qc_plots/scatter_plot_per_sample/",f,".pdf"), width=14)
	options(repr.plot.height = 6, repr.plot.width = 15)
	
	p1 <- all_samples[[f]]@meta.data %>% 
	  filter(percent.mouse < 5) %>% 
	  ggplot(aes(x = nCount_RNA, 
	             y = percent.mt)) +
	  ggtitle(f) + 
	  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
	  geom_hline(yintercept = 15, linetype = "dashed", color = "black", size = 1.5) + 
	  geom_vline(xintercept = 10000, linetype="dashed", color = "black", size = 1.5) +
	  theme_bw()

	p2 <-  all_samples[[f]]@meta.data %>% 
	  filter(percent.mouse < 5) %>% 
	  ggplot(aes(x = nCount_RNA, 
	             y = nFeature_RNA)) +
	  ggtitle(f) + 
	  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
		geom_hline(yintercept=2500, 
		           linetype="dashed", 
		           color = "black", 
		           size=1.5) + 
	  geom_vline(xintercept = 10000, 
	             linetype = "dashed", 
	             color = "black", 
	             size=1.5) +
	  theme_bw()

	print(p1+p2)

	### histograph
	options(repr.plot.height = 6, repr.plot.width = 24)
	plot1 <-  hist(all_samples[[f]]$nCount_RNA, breaks = 100, main=f, xlab = "Count depth (Counts/cell)", col='wheat')
	plot1 <-  abline(v=400,col="red")
	plot2 <-  hist(all_samples[[f]]$nFeature_RNA, breaks = 200, main=f, xlab = "Number of genes per cell", xlim = c(0, 4000), col='wheat')
	plot2 <-  abline(v=250,col="red")

	print(plot1 + plot2)

	dev.off()
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(paste0(output_dir,"/qc_plots/mito_ribo_stress_qc_final.pdf"), width=14, height=6)
for(s in names(all_samples)){
  mt_ribo_stress <- all_samples[[s]]@meta.data %>% filter(percent.mouse < 5)
  
  p1 <- ggplot(mt_ribo_stress, aes(x=percent.ribo, y=percent.mt))+geom_point(alpha=0.25, size=0.5)+ylim(0,100)+xlim(0,100)+ggtitle(s)
  p2 <- ggplot(mt_ribo_stress, aes(x=percent.stress, y=percent.mt))+geom_point(alpha=0.25, size=0.5)+ylim(0,100)+xlim(0,100)+ggtitle(s)
  
  # with marginal histogram
  p11 <- ggMarginal(p1, type="histogram")
  p21 <- ggMarginal(p2, type="histogram")
  print(as.ggplot(p11) + as.ggplot(p21))
  
}
dev.off()


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for(s in names(all_samples)){
    print(s)
    all_samples[[s]] <- RenameCells(all_samples[[s]], new.names=paste0(colnames(all_samples[[s]]),"_",s))
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
merged_samples <- merge(all_samples[[1]], all_samples[-1])
merged_samples

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p <- merged_samples@meta.data %>% 
  ggplot(aes(y = nFeature_RNA,
             x = percent.mouse)) + 
  geom_point(alpha = 0.1,
             size = 0.05) +
  geom_hline(yintercept = 200, linetype = 'dashed') +
  geom_vline(xintercept = c(5, 95), linetype = 'dashed') +
  scale_y_log10()

ggMarginal(p, type = 'densigram')


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
human_cells_to_keep <- merged_samples@meta.data %>% 
  filter(percent.mouse < 5 &
           nFeature_RNA > 200 & 
           nFeature_RNA < 6000 & 
           percent.mt > 1 &
           percent.mt < 20 & 
           percent.ribo > 0 &
           percent.ribo < 50 &
           percent.stress > 0 &
           percent.stress < 30) %>% 
  rownames()

mouse_cells_to_keep <- merged_samples@meta.data %>% 
  filter(percent.mouse > 95) %>% 
  rownames()

merged_samples <- subset(merged_samples, 
                         cells = c(human_cells_to_keep, mouse_cells_to_keep))

merged_samples

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
merged_samples@meta.data


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
adt <- merged_samples@meta.data %>% 
  select(id, nCount_ADT) %>% 
  group_by(id) %>% 
  summarize(median = median(nCount_ADT)) %>% 
  arrange(median) %>% 
  mutate(ADT = ifelse(median > 0, 'Yes', 'No')) %>% 
  select(id, ADT)

merged_samples@meta.data %>% left_join(adt) %>% 
  ggplot(aes(x = id,
             y = percent.stress)) +
  geom_violin() +
  facet_wrap(~ADT)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
merged_samples <- Seurat::NormalizeData(merged_samples, 
                                        assay = "RNA") %>% 
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 2000, 
                       assay = "RNA") %>%
  ScaleData(verbose = FALSE, 
            assay = "RNA") %>%
  RunPCA(npcs = 30, 
         verbose = FALSE, 
         assay = "RNA")


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
merged_samples %>% write_rds(file.path(data_dir, 'seuratobj_merged.rds'))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sessionInfo()

