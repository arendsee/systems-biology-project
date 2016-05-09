require(dplyr)
require(magrittr)
require(data.table)

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

#' Load a promoter file and remove rows with missing data
#' 
#' @export
#' @param p.path full filename for the promoter file
#' @return data.frame object containing promoter sites
#' @examples
#' d <- load_promoter('INPUT/athamap/HAHB4.txt')
load_promoter <- function(p.path){
  read.delim(p.path) %>% filter(!is.na(Relative.distance))
}

#' Add density and adjusted score columns to promoter table
#' 
#' @export
#' @param d data.frame object containing promoter sites 
#' @return d with new columns
#' @examples
#' d <- load_promoter('INPUT/athamap/HAHB4.txt')
#' d <- add_adjusted_scores(d)
add_adjusted_scores <- function(d){
  dens <- density(d$Relative.distance, n=1024)

  filter(d, Relative.orientation == '+') %>%
    mutate(density = interpolate(Relative.distance, dens$x, dens$y)) %>%
    mutate(adjusted_score = Score + 7 * log2(density / max(density)) ) %>%
    filter(adjusted_score > 0)
}
square_function <- function(d){
  filter(d, Relative.orientation == '+') %>%
    mutate(adjusted_score = ifelse(abs(Relative.distance) < 500, Score, 0)) %>%
    filter(adjusted_score > 0)
}

#' Group loci by model, collapsing on max adjusted_score
#' 
#' @export
#' @param d data.frame object containing promoter sites 
#' @return data.frame with two columns (locus, adjusted_score)
#' @examples
#' d <- load_promoter('INPUT/athamap/HAHB4.txt')
#' d <- add_adjusted_scores(d)
#' d <- reduce_to_locus(d)
reduce_to_locus <- function(d){
  mutate(d, locus = gsub('\\.\\d+', '', Gene)) %>%
    select(locus, adjusted_score) %>%
    group_by(locus) %>%
    summarise(adjusted_score=max(adjusted_score))
}

#' Build promoter matrix
#' 
#' @param promdir Directory containing all Athamap tables
#' @param method The 'square' method keeps only predictions ranging from -500 to 500.
#' @return Promoter matrix
#' @examples
#' m <- build_promoter_matrix('INPUT/athamap')
build_promoter_matrix <- function(promdir, method='square'){
  norms <- list()
  if(method == 'square'){
    adjust_score <- square_function
  } else {
    adjust_score <- add_adjusted_scores
  }
  for (f in list.files(promdir, pattern='*.txt', full.names=TRUE)){
    p.name <- sub('.txt', '', basename(f))
    d <- load_promoter(f)    %>%
         adjust_score %>%
         reduce_to_locus
    names(d)[2] <- p.name
    d <- data.table(d)
    setkey(d, locus)
    norms[[p.name]] <- d
  }
  d <- Reduce(function(...) merge(..., all=T, by='locus'), norms)
  d$locus <- gsub('At(\\d)g', 'AT\\1G', d$locus)
  d[is.na(d)] <- 0
  d
}
