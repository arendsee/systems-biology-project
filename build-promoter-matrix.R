#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(
  formatter_class='argparse.RawTextHelpFormatter',
  description='Build promoter score matrix based on Athamap tables. Writes matrix to STDOUT',
  usage='build-promoter-matrix.R <promoter_folder>')

parser$add_argument(
  'promdir',
  help='Folder containing all Athamap tables'
)

args <- parser$parse_args()

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("data.table"))

#' Interpolate densities for vector of values
#' 
#' @param a vector of values
#' @param x values for which a density has been calculated
#' @param y estimated density at position x
#' @return interpolated densities for vector a
#' @examples
#' d <- density(rnorm(1000), from=-3, to=3) 
#' interpolate(c(-1.234, 1.2322), d$x, d$y)
interpolate <- function(a, x, y){
  y[length(y)+1] <- y[length(y)] # just to avoid index errors latter
  step <- (max(x) - min(x)) / (length(x) - 1)
  i <- ceiling(pmax(0, a - min(x)) / step)
  y[i] + (y[i+1] - y[i]) * (a - x[i]) / step
}


#' Build promoter matrix
#' 
#' @param promdir Directory containing all Athamap tables
#' @return Promoter matrix
build_promoter_matrix <- function(promdir){
    norms <- list()
    for (f in list.files(promdir, pattern='*.txt')){
      promname <- sub('.txt', '', f)
      print(promname)

    d <- read.delim(file.path(promdir, f)) %>%
         filter(!is.na(Relative.distance))

    dens <- density(d$Relative.distance, n=1024)

    d <- filter(d, Relative.orientation == '+') %>%
         mutate(density = interpolate(Relative.distance, dens$x, dens$y)) %>%
         mutate(adjusted_score = Score + 7 * log2(density / max(dens$y))) %>%
         filter(adjusted_score > 0) %>%
         mutate(locus = gsub('\\.\\d+', '', Gene)) %>%
         select(locus, adjusted_score) %>%
         group_by(locus) %>%
         summarise(adjusted_score=max(adjusted_score))
    names(d)[2] <- promname

    d <- data.table(d)
    setkey(d, locus)

    norms[[promname]] <- d
  }

  d <- Reduce(function(...) merge(..., all=T, by='locus'), norms)
  d[is.na(d)] <- 0
  d <- as.data.frame(d)
  d$locus <- NULL
  rownames(d) <- gsub('At(\\d)g', 'AT\\1G', rownames(d))
  as.matrix(d)
}


d <- build_promoter_matrix(args$promdir)

write.table(as.matrix(d))


# d <- read.table("data/tf-counts.tab")
# names(d) <- c('model', 'tf', 'count')
#
# adj2net <- function(adjmatrix, ...){
#     graph_from_adjacency_matrix(adjmatrix, mode="undirected", diag=FALSE, ...)
# }
# # Return only the largest connected component
# largest_component <- function(x){
#     components(x) %$%
#         which(membership != which.max(csize)) %>%
#         delete_vertices(graph=x)
# }
# # Remove components of size less than or equal to k
# prune <- function(g, k=1){
#     components(g) %$%
#         which(membership %in% which(csize <= k)) %>%
#         delete_vertices(graph=g)
# }
#
# nrow(d)
# nlevels(d$tf)
#
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
