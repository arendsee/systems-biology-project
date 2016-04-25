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

parser$add_argument(
  'output',
  help='Output matrix filename'
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
  d$locus <- gsub('At(\\d)g', 'AT\\1G', d$locus)
  d
}

d <- build_promoter_matrix(args$promdir)

write.table(d, file=args$output, row.names=FALSE, quote=FALSE, sep=" ", col.names=TRUE)
