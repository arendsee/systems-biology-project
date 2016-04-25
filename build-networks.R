#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(
  formatter_class='argparse.RawTextHelpFormatter',
  description='Do cool stuff with an expression and promoter matrix',
  usage='build-networks.R')

parser$add_argument(
  '-v', '--version',
  action='store_true',
  default=FALSE)

parser$add_argument(
  '-e', '--expression_matrix',
  help='Expression matrix (output of build-expression-matrix.R)'
)

parser$add_argument(
  '-p', '--promoter_matrix',
  help='Promoter similarity matrix (output of build-promoter-matrix.R)'
)

args <- parser$parse_args()

require(tidyr)
require(reshape2)
require(igraph)
require(dplyr)
require(minet)

  # TODO - remove this argument bypass
  args <- list(
    expression_matrix='matrices/expression.mat',
    promoter_matrix='matrices/promoter.mat'
  )

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

em <- read.table(args$expression_matrix, header=TRUE, row.names=1)
pm <- read.table(args$promoter_matrix, header=TRUE, row.names=1)
orphans <- read.table('data/orphans-locus', stringsAsFactors=F)$V1

loci <- intersect(rownames(em), rownames(pm))

incl_orphans <- intersect(loci, orphans)
em <- em[loci, ]
pm <- pm[loci, ]

correlation.between.promoters <- function(){
  mim <- build.mim(pm, estimator="spearman")
  aracne.net <- aracne(mim)
  ag <- adj2net(aracne.net > 0.3)
  plot(ag, vertex.size=1)
}

correlation.between.conditions <- function(){
  mim <- build.mim(em, estimator="spearman")
  aracne.net <- aracne(mim)
  ag <- adj2net(aracne.net > 0.2)
  plot(ag, vertex.size=1)
}

#' Build orphan promoter network
#' 
#' @param k keep only motifs with less than n/k non-zero scores
pm.net <- function(k=5, cor.cutoff=0.55, cor.method='spearman'){

  spm <- pm[incl_orphans, ]
  spm <- spm[, which(colSums(spm > 0) < (nrow(spm) / k))]
  spm[, 'TaNAC69.1.'] <- NULL
  spm <- spm[which(rowSums(spm) != 0), ]
  spm <- sapply(spm, function(x) x / max(x))
  p.cor <- cor(t(spm), method=cor.method)
  ag <- adj2net(p.cor > cor.cutoff)
  plot(ag, vertex.label=NA, vertex.size=0.2)

}

em.net <- function(){

  sem <- em[incl_orphans, ]
  sem <- sem[, -c(7,8,9,10)] # these ones are WAY too correlated
  sem[is.na(sem)] <- 0
  sem <- spm[which(rowSums(spm) != 0), ]
  mim <- build.mim(t(sem), estimator="spearman")
  aracne.net <- aracne(mim)
  ag <- adj2net(aracne.net > 0.4)
  plot(ag, vertex.label=NA, vertex.size=1)

  p.cor <- cor(t(sem), method='spearman')
  ag <- adj2net(p.cor > 0.6)
  plot(ag, vertex.label=NA, vertex.size=0.2)
}

# plot(ag, vertex.label=NA, vertex.size=1)

# # aracne_eps <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
# # ag <- lapply(aracne_eps, function(i) adj2net(aracne.net > i))
# # par(mfrow=c(3,2), mar=c(0,0,0,0))
# # lapply(ag, largest_component) %>%
# #     lapply(plot, vertex.label=NA, vertex.size=1)
#
# ag <- adj2net(aracne.net > 0.9)
# plot(ag, vertex.label=NA, vertex.size=1)


# hist(nd$freq, breaks=20)
#
#   nd <- as.character(d$tf) %>% count
#   uncommon <- subset(nd, freq < 10000) %$% as.character(x)
#
#   d <- d[d$tf %in% uncommon, ] %>% droplevels
#
#   nrow(d)
#   nlevels(d$tf)
#
# m  <- acast(d, model ~ tf, sum, fill=0)
# tf.cor <- cor(m, method="spearson")
# par(mfrow=c(1,1), mar=c(0,0,0,0))
# adj2net(tf.cor > 0.7) %>%
#   # plot(vertex.size=1)
#   plot(vertex.label=NA, vertex.size=1)
#
# m.sam <- m[sample.int(nrow(m), 10000), ]
#
# gen.cor <- cor(t(m.sam), method="spearman")
# par(mfrow=c(1,1), mar=c(0,0,0,0))
# adj2net(gen.cor > 0.8) %>%
#   prune(k=2) %>%
#   plot(vertex.label=NA, vertex.size=1)
#   # plot(vertex.size=1)
#
# # m <- as.matrix(d[c(1,2)]) %>%
# #   graph_from_edgelist(directed=FALSE) %>%
# #   get.adjacency
#
# # g <- make_ego_graph(g, 1, nodes=V(g)[V(g)$name %in% uncommon]) %>%
# #   lapply(function(x) connect.neighborhood(x, 2))
#
# # a <- data.frame(x=letters[sample.int(3,10, rep=T)],
# #                 y=letters[sample.int(3,10, rep=T)])
