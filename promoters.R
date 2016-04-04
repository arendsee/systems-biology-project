require(plyr)
require(magrittr)
require(igraph)
require(reshape2)
require(tidyr)
require(bigmemory)
require(bigalgebra)
require(bigpca)


d <- read.table("tf-counts.tab")
names(d) <- c('model', 'tf', 'count')

nrow(d)
nlevels(d$tf)

nd <- as.character(d$tf) %>% count
hist(nd$freq, breaks=20)

  uncommon <- subset(nd, freq < 10000) %$% as.character(x)

  d <- d[d$tf %in% uncommon, ] %>% droplevels

  nrow(d)
  nlevels(d$tf)

m  <- acast(d, model ~ tf, sum, fill=0) > 0
mt <- t(m)
mt <- as.big.matrix(mt)
m  <- as.big.matrix(m)
b  <- m %*% mt

# m <- as.matrix(d[c(1,2)]) %>%
#   graph_from_edgelist(directed=FALSE) %>%
#   get.adjacency

# g <- make_ego_graph(g, 1, nodes=V(g)[V(g)$name %in% uncommon]) %>%
#   lapply(function(x) connect.neighborhood(x, 2))

# a <- data.frame(x=letters[sample.int(3,10, rep=T)],
#                 y=letters[sample.int(3,10, rep=T)])
