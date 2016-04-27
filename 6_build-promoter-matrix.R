#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(
  formatter_class='argparse.RawTextHelpFormatter',
  description='Build promoter score matrix based on Athamap tables. Writes matrix to STDOUT',
  usage='build-promoter-matrix.R <promoter_folder>')

parser$add_argument(
  '-o', '--output',
  help='Output matrix filename',
  default='OUTPUT/matrices/promoter.mat'
)

parser$add_argument(
  '-p', '--promdir',
  help='Folder containing all Athamap tables',
  default='INPUT/athamap'
)

args <- parser$parse_args()

source('promoter-normalization.R')

d <- build_promoter_matrix(args$promdir)

write.table(d, file=args$output, row.names=FALSE, quote=FALSE, sep=" ", col.names=TRUE)
