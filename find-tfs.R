
# This script is designed to be run on a cluster with lots of memory (>16G).
# It requires the following files be present

output <- 'tf-orf-output'

require(igraph)
require(minet)

adj2net <- function(adjmatrix, ...){
    graph_from_adjacency_matrix(adjmatrix, mode="undirected", diag=FALSE, ...)
}

# load expression data
em <- read.table('expression.mat', header=TRUE, row.names=1)
orphans <- read.table('orphans', stringsAsFactors=F)$V1

tf.info <- read.delim('tf.tab')

orphans <- intersect(orphans, rownames(em))
tfs     <- intersect(tf.info$locus, rownames(em))

loci <- Reduce(union, list(rownames(em), orphans , tf.info$locus))

    # DELETE THIS SECTION
    long.tf.list <- read.delim('long-tf-list.tab')$gene_model %>%
        gsub(pattern='\\.\\d+', replacement='') %>% unique
    loci <- Reduce(union, list(orphans, tf.info$locus, long.tf.list, sample(rownames(em), 2000)))

em <- em[loci, ]
em[is.na(em)] <- 0

# build ARACNE graph
mim <- build.mim(t(em), estimator="spearman")
aracne.net <- aracne(mim)

others  <- setdiff(loci, c(orphans, tfs, long.tf.list))

testset <- list(
    orphans = orphans,
    r1      = sample(others, length(orphans)),
    r2      = sample(others, length(orphans)),
    r3      = sample(others, length(orphans)),
    r4      = sample(others, length(orphans)),
    r5      = sample(others, length(orphans))
)


cutoffs <- seq(from=0.5, to=0.975, by=0.025)
g.degrees  = list()
g.stats <- data.frame(
    cutoff   = cutoffs,
    avg_path = rep(0, length(cutoffs)),
    density  = rep(0, length(cutoffs))
)
for (i in 1:length(cutoffs)){
  dname <- file.path(output, paste0('c', cutoffs[i]))
  dir.create(dname, recursive=TRUE, showWarnings=FALSE)
  ag <- adj2net(aracne.net > cutoffs[i])
  g.stats[i, 'avg_path'] <- mean_distance(ag)
  g.stats[i, 'density'] <- graph.density(ag)
  g.degrees[[i ]] <- data.frame(degree=degree(ag)) 
  g.degrees[[i ]]$cutoff <- cutoffs[i]
  for (n in names(testset)){
      which(!vertex_attr(ag)$name %in% c(tfs, testset[[n]])) %>%
          delete_vertices(graph=ag) %>%
          write.graph(file.path(dname, paste0(n, '.mat')), 'ncol')
  }
}

g.degrees <- Reduce(rbind, g.degrees)
write.table(g.degrees, file.path(output, 'graph-degrees.tab'), sep="\t", row.names=FALSE)

write.table(g.stats, file.path(output, 'graph-stats.tab'), sep="\t", row.names=FALSE)
