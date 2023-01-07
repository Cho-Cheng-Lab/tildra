# integrating single-cell spatial RNA-seq data

# (!) ref1: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_07_spatial.html#Identification_of_Spatially_Variable_Features
# (!) ref2: https://satijalab.org/seurat/articles/spatial_vignette.html#working-with-multiple-slices-in-seurat-1

# load the libraries
library(Seurat) 
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

# specify input, output and metadata
output_dir <- "analysis/output/10_visium/"
samples_dir <- "analysis/visium_samples/"
metadata_df <- read.table("/media/pk3/143F2A3651271F75/special_projects/1_scrnaseq_ucsf_jeffrey_analysis_1/BI/tildra/analysis/output/00_process/metadata_subset.tsv", header=TRUE, sep='\t', stringsAsFactors=FALSE)

setwd(samples_dir)

all_samples <- list()

# I renamed "tissue_positions.csv" to "tissue_positions_list.csv"
# load in the samples and remove genes with all zero counts: https://github.com/satijalab/sctransform/issues/86
for(s in dir(samples_dir)){
	print(s)
	all_samples[[s]] <- Load10X_Spatial(data.dir=paste0(samples_dir,s), filename="filtered_feature_bc_matrix.h5", slice=metadata_df[metadata_df$id == s,]$id)
	all_samples[[s]] <- subset(all_samples[[s]], subset = nCount_Spatial > 0)
}

# add metadata
for(n in names(all_samples)){
	md <- subset(metadata_df, metadata_df$id == n)
	
	for(m in colnames(md)[c(1:5,9)]){
		print(paste0("adding ",m))
		all_samples[[n]][[m]] <- md[[m]]
	}
}

# run SCT on all datasets
all_samples = lapply(all_samples, SCTransform, assay = "Spatial", method = "poisson")

# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB

# select integration features and run SCTIntegration
features = SelectIntegrationFeatures(all_samples, nfeatures = 3000, verbose = FALSE)
all_samples <- PrepSCTIntegration(object.list = all_samples, anchor.features = features,
    verbose = FALSE)

# samples integration
int.anchors <- FindIntegrationAnchors(object.list = all_samples, normalization.method = "SCT",
    verbose = FALSE, anchor.features = features)
samples_integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
    verbose = FALSE)

# dimensionality reduction and clustering
samples_integrated <- RunPCA(samples_integrated, verbose = FALSE)
samples_integrated <- FindNeighbors(samples_integrated, dims = 1:30)
samples_integrated <- FindClusters(samples_integrated, resolution = seq(0.1,2,by=0.1), verbose = FALSE)
samples_integrated <- RunUMAP(samples_integrated, dims = 1:30)

# (!) https://stackoverflow.com/questions/73131436/error-in-funleft-right-non-numeric-argument-to-binary-operator-when-runni
for(n in names(samples_integrated@images)){
	samples_integrated@images[[n]]@coordinates[["tissue"]] <- as.integer(samples_integrated@images[[n]]@coordinates[["tissue"]])
	samples_integrated@images[[n]]@coordinates[["row"]] <- as.integer(samples_integrated@images[[n]]@coordinates[["row"]])
	samples_integrated@images[[n]]@coordinates[["col"]] <- as.integer(samples_integrated@images[[n]]@coordinates[["col"]])
	samples_integrated@images[[n]]@coordinates[["imagerow"]] <- as.integer(samples_integrated@images[[n]]@coordinates[["imagerow"]])
	samples_integrated@images[[n]]@coordinates[["imagecol"]] <- as.integer(samples_integrated@images[[n]]@coordinates[["imagecol"]])
}

# (!) to choose one image (or slice): https://github.com/satijalab/seurat/issues/5774
# p1 <- SpatialDimPlot(samples_integrated, images = c("Skin234","Skin244"))+NoLegend()
# p2 <- SpatialDimPlot(samples_integrated, images = c("skin296","skin269"))+NoLegend()
# p3 <- SpatialDimPlot(samples_integrated, images = c("Skin235","skin279"))+NoLegend()
# p4 <- SpatialDimPlot(samples_integrated, images = c("skin285","skin292"))+NoLegend()
# p5 <- SpatialDimPlot(samples_integrated, images = c("skin289","skin290"))+NoLegend()
# p6 <- SpatialDimPlot(samples_integrated, images = c("skin305","skin304"))+NoLegend()

# p1 / p2 / p3 / p4 / p5 / p6

Idents(samples_integrated) <- "integrated_snn_res.1"

pdf(paste0("analysis/output/10_visium/UMAP_integrated_per_cluster_and_treatment.pdf"), height=8, width=14)
DimPlot(samples_integrated, reduction = "umap", group.by = c("ident", "treatment"))
dev.off()

for(n in unique(samples_integrated$id)){
	pdf(paste0("analysis/output/10_visium/",n,"_spatial_umap.pdf"), height=10, width=10)
	print(SpatialDimPlot(samples_integrated, images = n, alpha=1.5))
	dev.off()
}

# DE between pre vs post - resolution 1
clusters <- as.character(unique(Idents(samples_integrated)))

OUT <- createWorkbook()
for(cluster in clusters){
	tryCatch({
	print(paste0(cluster))
	markers <- FindMarkers(samples_integrated, ident.1=cluster, logfc.threshold=0, 
	min.pct=0.01, test.use="MAST") %>% 
	mutate(., gene = rownames(.))

	addWorksheet(OUT, paste0(cluster))
	writeData(OUT, sheet = paste0(cluster), x = markers)
	}, error=function(e){message("Cell group 1 has fewer than 3 cells")})
}

saveWorkbook(OUT, paste0("analysis/output/11_visium_de/marker_genes_per_cluster.xlsx"))

# Annotation based on the markers above - marker_genes_per_cluster.xlsx
new.cluster.ids <- c("Fibroblast","Keratinocytes","Fibroblast","Endothelial","Fibroblast","Fibroblast",
	"Fibroblast","Fibroblast","Fibroblast","Keratinocytes","Fibroblast","Fibroblast",
	"Keratinocytes","Keratinocytes","Fibroblast","Fibroblast","Keratinocytes","Keratinocytes")

samples_integrated$annotation <- Idents(samples_integrated)

names(new.cluster.ids) <- levels(samples_integrated)
samples_integrated <- RenameIdents(samples_integrated, new.cluster.ids)

samples_integrated$anno_res_1 <- Idents(samples_integrated)

# UMAP with annotation
pdf("analysis/output/10_visium/UMAP_with_annotation.pdf",width=8, height=8)
DimPlot(samples_integrated)
dev.off()

# save the object
#saveRDS(samples_integrated, "/media/pk3/143F2A3651271F75/special_projects/8_single_cell_spatial_jeffrey/integrated_samples.rds")


# extract data - as per David's request
## get the donor, treatment and cluster columns for all cells
table1 <- samples_integrated@meta.data[,c("id","treatment","integrated_snn_res.1","anno_res_1")]
table1 <- data.frame(cell_barcode = rownames(table1), donor=table1$id, treatment=table1$treatment, cluster=table1$integrated_snn_res.1, annotation=table1$anno_res_1)
write.table(table1, "analysis/output/10_visium/cell_barcodes_donor_treatment_and_cluster_at_res_1.txt", 
	col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')

## get the normalized counts for the genes of interest
normalized_counts_specific_genes <- as.matrix(samples_integrated@assays$SCT@data) %>% subset(., rownames(.) %in% c("IL17A","IL17F","IL23A","IL23R","ZFP36L2","ZFP36"))
normalized_counts_specific_genes <- cbind(genes=rownames(normalized_counts_specific_genes), normalized_counts_specific_genes)
rownames(normalized_counts_specific_genes) <- NULL
normalized_counts_specific_genes <- as.data.frame(normalized_counts_specific_genes)

write.table(normalized_counts_specific_genes, "analysis/output/10_visium/normalized_counts_for_genes_of_interest.txt", 
	col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')
