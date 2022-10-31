# single cell spatial analysis with 'denoised' normalized counts from SPROD
library(Seurat)
library(Matrix)
library(metap) 
library(MAST) 
library(scDblFinder) 
library(ggplot2) 
library(cowplot) 
library(patchwork) 
library(writexl)
library(dplyr)
library(xlsx)
library(harmony)
library(openxlsx)
library(rgl)
library(ggExtra)
library(ggplotify)


output_dir <- "tildra/output/10_visium/"
samples_dir <- "tildra/analysis/input/visium_samples/"
metadata_df <- read.table("tildra/analysis/output/00_process/metadata_subset.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)

setwd(samples_dir)

all_samples <- list()

# I renamed "tissue_positions.csv" to "tissue_positions_list.csv"
# load in the samples and remove genes with all zero counts: https://github.com/satijalab/sctransform/issues/86
for(s in dir(samples_dir)){
	print(s)
	all_samples[[s]] <- Load10X_Spatial(data.dir=paste0(samples_dir,s), filename="filtered_feature_bc_matrix.h5", slice=metadata_df[metadata_df$donor == s,]$donor)
	all_samples[[s]] <- subset(all_samples[[s]], subset = nCount_Spatial > 0)
}

# add metadata
for(n in names(all_samples)){
	md <- subset(metadata_df, metadata_df$donor == n)
	
	for(m in colnames(md)){
		print(paste0("adding ",m))
		all_samples[[n]][[m]] <- md[[m]]
	}
}

all_samples = lapply(all_samples, SCTransform, assay = "Spatial", method = "poisson")


# (!) https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
for(n in names(all_samples)){
	all_samples[[n]]@images[[n]]@coordinates[["tissue"]] <- as.integer(all_samples[[n]]@images[[n]]@coordinates[["tissue"]])
	all_samples[[n]]@images[[n]]@coordinates[["row"]] <- as.integer(all_samples[[n]]@images[[n]]@coordinates[["row"]])
	all_samples[[n]]@images[[n]]@coordinates[["col"]] <- as.integer(all_samples[[n]]@images[[n]]@coordinates[["col"]])
	all_samples[[n]]@images[[n]]@coordinates[["imagerow"]] <- as.integer(all_samples[[n]]@images[[n]]@coordinates[["imagerow"]])
	all_samples[[n]]@images[[n]]@coordinates[["imagecol"]] <- as.integer(all_samples[[n]]@images[[n]]@coordinates[["imagecol"]])
}

# (!) Use this block of code fro line 59 to 81 if I am running the pipeline from the beginning
all_samples_denoised <- all_samples

pdf(paste0(output_dir,"/KRT10_and_COL1A1_spatial_expression_before_and_after_denoising.pdf"), width=12, height=6)
for(n in names(all_samples)){
	cat(n,"\n")
	p1 <- SpatialFeaturePlot(all_samples[[n]], features=c("KRT10","COL1A1"))+ggtitle(paste0(n," Before denoising"))

	# here is the slot of the normalized counts: skin235@assays$SCT@data
	## read in the denoised matrix
	my_matrix_zero <- t(read.table(paste0(samples_dir,n,"/sprod_Denoised_matrix.txt"), 
		header=TRUE, sep='\t', stringsAsFactors=FALSE))

	# transform it to a sparse matrix
	my_matrix <- as.sparse(my_matrix_zero)

	all_samples_denoised[[n]]@assays$SCT@data <- my_matrix

	p2 <- SpatialFeaturePlot(all_samples_denoised[[n]], features=c("KRT10","COL1A1"))+ggtitle(paste0(n,"After denoising"))

	print(p1 | p2)
}

dev.off()

# saving the objects all_samples and all_samples_denoised
# save(all_samples, all_samples_denoised, file="/media/kurdiabed/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/all_samples_and_all_samples_denoised_objects.RData")

# (!) use this block of code just if I need to re-plot, after loading all_samples and all_samples_denoised objects
# just load "all_samples_and_all_samples_denoised_objects.RData" file and do the plotting

pdf(paste0(output_dir,"/KRT10_and_COL1A1_spatial_expression_before_and_after_denoising.pdf"), width=12, height=6)
for(n in names(all_samples)){
	cat(n,"\n")
	p1 <- SpatialFeaturePlot(all_samples[[n]], features=c("KRT10","COL1A1"), image.alpha=0)+ggtitle(paste0(n," Before denoising"))
	p2 <- SpatialFeaturePlot(all_samples_denoised[[n]], features=c("KRT10","COL1A1"), image.alpha=0)+ggtitle(paste0(n,"After denoising"))
	print(p1 | p2)
}
dev.off()
