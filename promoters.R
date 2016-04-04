require(plyr)
require(magrittr)
require(igraph)
require(reshape2)
require(tidyr)

d <- read.table("data/tf-counts.tab")
names(d) <- c('model', 'tf', 'count')

adj2net <- function(adjmatrix, ...){
    graph_from_adjacency_matrix(adjmatrix, mode="undirected", diag=FALSE, ...)
}
# Return only the largest connected component
largest_component <- function(x){
    components(x) %$%
        which(membership != which.max(csize)) %>%
        delete_vertices(graph=x)
}
# Remove components of size less than or equal to k
prune <- function(g, k=1){
    components(g) %$%
        which(membership %in% which(csize <= k)) %>%
        delete_vertices(graph=g)
}

nrow(d)
nlevels(d$tf)

hist(nd$freq, breaks=20)

  nd <- as.character(d$tf) %>% count
  uncommon <- subset(nd, freq < 10000) %$% as.character(x)

  d <- d[d$tf %in% uncommon, ] %>% droplevels

  nrow(d)
  nlevels(d$tf)

m  <- acast(d, model ~ tf, sum, fill=0)
tf.cor <- cor(m, method="spearson")
par(mfrow=c(1,1), mar=c(0,0,0,0))
adj2net(tf.cor > 0.7) %>%
  # plot(vertex.size=1)
  plot(vertex.label=NA, vertex.size=1)

m.sam <- m[sample.int(nrow(m), 10000), ]

gen.cor <- cor(t(m.sam), method="spearman")
par(mfrow=c(1,1), mar=c(0,0,0,0))
adj2net(gen.cor > 0.8) %>%
  prune(k=2) %>%
  plot(vertex.label=NA, vertex.size=1)
  # plot(vertex.size=1)

# m <- as.matrix(d[c(1,2)]) %>%
#   graph_from_edgelist(directed=FALSE) %>%
#   get.adjacency

# g <- make_ego_graph(g, 1, nodes=V(g)[V(g)$name %in% uncommon]) %>%
#   lapply(function(x) connect.neighborhood(x, 2))

# a <- data.frame(x=letters[sample.int(3,10, rep=T)],
#                 y=letters[sample.int(3,10, rep=T)])
