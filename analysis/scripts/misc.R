```{r}
library(ape)
seuratobj <- BuildClusterTree(seuratobj, reduction = 'umap', dims = 1:2, assay = 'RNA')
tree <- as.hclust(Tool(seuratobj, slot = 'BuildClusterTree'))
```
```{r}
tree_height <- 100
plot(tree)
abline(h = tree_height, col = 'red')
```
### Merged clusters from phylogenetic tree
```{r}
cuts <- cutree(tree, k = 7)

tree_mapping <- tibble('seurat_clusters' = names(cuts),
                       'broad_clusters' = factor(cuts))

broad_clusters <- seuratobj@meta.data %>% rownames_to_column() %>% select(seurat_clusters) %>% left_join(tree_mapping)

seuratobj$broad_clusters <- broad_clusters$broad_clusters

p <- seurat_feature(seuratobj, feature = 'broad_clusters', facet_hide = TRUE, rasterize_dpi = 600, title = 'Broad Clusters')

ggsave(plot = p,
       filename = 'broad_clusters.png',
       h = 5,
       w = 5,
       dpi = 600,
       path = output_dir)

p
```