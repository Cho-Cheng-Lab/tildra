dims <- 1:20
neighbors <- 30
metric <- 'cosine'
prune <- 1/15
seed <- 123

set.seed(123)

test <- seuratobj %>% 
  FindNeighbors(reduction = 'harmony', 
                dims = dims, 
                k.param = neighbors, 
                annoy.metric = metric, 
                prune.SNN = prune, 
                seed = seed)

test@reductions$pca@cell.embeddings[1:5, 1:5]

test@graphs$RNA_nn %>% sum()
test@graphs$RNA_snn %>% sum()
