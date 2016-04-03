# source("https://bioconductor.org/biocLite.R")
# biocLite("SRAdb")

library(SRAdb)
library(magrittr)

sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

sra_con <- dbConnect(SQLite(),sqlfile)

tables <- function() dbListTables(sra_con)
fields <- function(table) dbListFields(sra_con, table)
query  <- function(cmd) dbGetQuery(src_con, cmd)

fetch <- function(fields, table, condition, limit){
  fields = paste(fields, collapse=", ")
  cmd=sprintf("select %s from %s", fields, table)
  if(!missing(condition)){
    condition = paste(condition, collapse=" AND ")
    cmd=sprintf("%s where %s", cmd, condition)
  }
  if(!missing(limit)){
    cmd=sprintf("%s limit %d", cmd, limit) 
  }
  dbGetQuery(sra_con, cmd)
}

# taxon_id:3702 library_strategy:RNA library_source:TRANSCRIPTOMIC library_layout:PAIRED

fetch("abstract", "study", "taxon_id = 3702", 5)

# fetch(
#   fields    = "study_title",
#   table     = "sra_ft", 
#   limit     = 5,
#   condition = c(
#     "taxon_id = '3702'",
#     "library_strategy = 'RNA-Seq'",
#     "library_source = 'TRANSCRIPTOMIC'",
#     "library_layout = 'PAIRED'"
#   )
# )



# Take all TAIR10 models

# Merge genes with identical coding sequences

# Build kallisto index with these genes

# Select set of interesting Arabidopsis studies
# - want wide variety of conditions
# - use SRAdb to get list of runids

# For each, fetch | fastq-dump | kallisto

# Merge all into final expression table
