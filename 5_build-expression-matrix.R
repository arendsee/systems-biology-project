#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(
  formatter_class='argparse.RawTextHelpFormatter',
  description='Build an expression matrix from Kallisto output',
  usage='build-expression-matrix <kallisto-output-directory>'
)

parser$add_argument(
  'output',
  help='Output matrix file'
)

parser$add_argument(
  '-k', '--kallisto-out',
  help='Directory of Kallisto output files (folders with SRA runid names)',
  default='INTERMEDIATE/kallisto-runs'
)

parser$add_argument(
  '-c' , '--conditions',
  help='Experimental design file with fields "run_accession" and  "condition"',
  default='INPUT/conditions.tab'
)

args <- parser$parse_args()


# # Install rhdf5 with the following code, if needed
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
#
# # Install sleuth with the following if needed
# source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")    # only if devtools not yet installed
# biocLite("pachterlab/sleuth")

require(sleuth)
require(magrittr)
require(rhdf5)
require(dplyr)
require(data.table)


#' Load experimental design
#' 
#' @export
#' @param data_dir
#' @param condition_file
#' @param sample_colname
#' @param condition_colname
#' @return desc
#' @examples
#' 
get_full_meta_study_ <- function(
    data_dir=args$kallisto_out,
    condition_file=args$conditions,
    sample_colname='run_accession',
    condition_colname='condition'
){
  base_dir=file.path(getwd(), data_dir)
  d <- read.delim(condition_file)
  d <- d[, c(sample_colname, condition_colname)] 
  names(d) <- c('sample', 'condition')
  d$path <- sapply(d$sample, function(id) file.path(base_dir, id))
  d
}



get_genemap_ <- function(d){
  transcripts <- read.delim(file.path(d$path[1], 'abundance.tsv'), stringsAsFactors=F)[[1]]
  data.frame(
    target_id=transcripts,
    locus=sub('\\.\\d+', '', transcripts)
  )
}



#' Run all Kallisto data through Sleuth
#' 
#' @return results table
process_sleuth_results_ <- function(full_data, control, treatment, genemap){
  d <- subset(full_data, condition %in% c(control, treatment)) %>% droplevels
  levels(d$condition) <- c(treatment, control)

  condition <- paste0('condition', control)

  d <- sleuth_prep(d, ~ condition, target_mapping=genemap) %>%
    sleuth_fit %>%
    sleuth_wt(which_beta=condition) %>%
    sleuth_results(condition)
}

fulldat <- get_full_meta_study_()
genemap <- get_genemap_(fulldat)


#' A wrapper for process_sleuth_results_
sleu <- function(...){
  process_sleuth_results_(full_data=fulldat, genemap=genemap, ...)
}

studies <- list(
   copper_shoot = sleu(control='D-shoot-Copper-deficient', treatment='shoot-Copper-deficient'),
   copper_root  = sleu(control='D-root-Copper-deficient',  treatment='root-Copper-deficient'),
   nitrate      = sleu(control='D-root_high_nitrate',      treatment='root_high_nitrate'),
   FM           = sleu(control='IM',                       treatment='FM'),
   ST3          = sleu(control='IM',                       treatment='ST3'),
   GA           = sleu(control='D-GA-BR',                  treatment='low-GA'),
   BR           = sleu(control='D-GA-BR',                  treatment='low-BR'),
   GA_BR        = sleu(control='D-GA-BR',                  treatment='low-GA_low-BR'),
   mold_1       = sleu(control='Waco9_0ddpi',              treatment='Waco9_1ddpi'),
   mold_3       = sleu(control='Waco9_0ddpi',              treatment='Waco9_3ddpi'),
   mold_5       = sleu(control='Waco9_0ddpi',              treatment='Waco9_5ddpi'),
   seed_FA_1    = sleu(control='fae1_7-8',                 treatment='fae1_9-10'),
   seed_FA_2    = sleu(control='fae1_7-8',                 treatment='fae1_11-12'),
   iron         = sleu(control='D-iron-deficiency',        treatment='iron-deficiency'),
   nematode_4   = sleu(control='D-nematode-4dpi',          treatment='nematode-4dpi'),
   nematode_10  = sleu(control='D-nematode-10dpi',         treatment='nematode-10dpi')
)
    
clean_to_fold <- function(a){
  group_by(a, locus) %>%
    slice(which.min(pval)) %>%
    ungroup %>%
    select(locus, b) %>%
    as.data.table %>%
    setkey(locus)
}

d <- lapply(studies, clean_to_fold)
d <- lapply(names(d), function(n) {names(d[[n]]) <- c('locus', n); d[[n]]})
d <- Reduce(function(...) merge(..., all=T, by='locus'), d)
d[is.na(d)] <- 0

write.table(d, file=args$output, row.names=FALSE, quote=FALSE, sep=" ", col.names=TRUE)
