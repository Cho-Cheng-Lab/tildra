mouse_genes <- rownames(test_mat) %>% str_detect('mm10')

mouse_mat <- test_mat[mouse_genes,]

mouse_cells <- colSums(mouse_mat) > 1000 

tibble(mouse_cells) %>% 
  ggplot(aes(x = mouse_cells + 0.001)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  geom_vline(xintercept = 1000)

mouse_mat_subset <- mouse_mat[,mouse_cells]

mouse_so <- CreateSeuratObject(counts = mouse_mat_subset)
mouse_so <- PercentageFeatureSet(mouse_so, pattern = 'Rp[sl]|Mrp[sl]', col.name = 'pct.ribo')

gene_sums <- rowSums(mouse_mat_subset)

gene_table <- tibble('gene' = names(gene_sums),
                     'count' = gene_sums) %>% 
  mutate(type = case_when(
    str_detect(gene, 'Rp[sl]') ~ 'ribo',
    str_detect(gene, 'Mrp[sl]') ~ 'ribo',
    str_detect(gene, 'mt-') ~ 'mito')) %>% 
  arrange(-count) %>% 
  mutate(rank = row_number()) %>% 
  drop_na() %>% 
  group_by(type) %>% 
  mutate(cumulative = cumsum(count),
         cumulative_f = cumulative/sum(count))
gene_table
