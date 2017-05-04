### Creating plots for all networks
devtools::install_github('thomasp85/ggforce')
devtools::install_github('thomasp85/ggraph')

library(ggraph)
library(igraph)
library(cowplot)

nets.sym <- lapply(paste0("./nets/symtrans/net_", 1:1000, "_symtrans.csv"), read.csv)
nets.asym <- lapply(paste0("./nets/asymtrans/net_", 1:1000, "_asymtrans.csv"), read.csv)

## Testing ggraph
teste <- graph.adjacency(as.matrix(nets.sym[[5]]), mode = "directed", diag = FALSE)
V(teste)$degree <- degree(teste, mode = 'in')

ggraph(teste, 'igraph', algorithm = 'kk') +
    geom_edge_fan() +
    geom_node_point(aes(size = degree), color = "red")
