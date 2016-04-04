# source("https://bioconductor.org/biocLite.R")
# biocLite("SRAdb")

library(SRAdb)
library(magrittr)

sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists('SRAmetadb.sqlite')) sqlfile <<- getSRAdbFile()

sra_con <- dbConnect(SQLite(),sqlfile)

tables <- function() dbListTables(sra_con)
fields <- function(table) dbListFields(sra_con, table)
query  <- function(cmd) dbGetQuery(conn=sra_con, statement=cmd)

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
  dbGetQuery(conn=sra_con, statement=cmd)
}

stuid = c(
  'SRP009850', # shoot apical meristems stages
  'SRP018034', # leaf senesence stages
  'SRP053394', # hi/lo temp
  'SRP059384', # roots versus shoots
  'SRP059724', # (hi/lo Pi) X (mock/Ct/Ci) X (6,10,16,24 dpi)
  'SRP063314', # leaf, flower, root
  'SRP064782', # circadian clock ???
  'SRP063421', # (hi/lo Pi) X (root, shoot) X (long, short exposure)
  'SRP069266'  # high altitude adaptation
)

stuid2 = c(
  'SRP033660', # (17 accession) X (seedling, root, flower) X (3 bio reps)
  'SRP036643'  # (160 Swedish accession) X (10C, 16C)
)

fjoin <- function(x) sprintf("'%s'", paste0(x, collapse="', '"))

getRuns <- function(s){
  sprintf("'%s'", paste0(s, collapse="', '")) %>%
    sprintf(fmt='select * from sra_ft where study_accession in ( %s )') %>%
    query()
}

runs  <- getRuns(stuid)[c()]
runs2 <- getRuns(stuid2)


# taxon_id:3702 library_strategy:RNA library_source:TRANSCRIPTOMIC library_layout:PAIRED


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

# ---------

# In Sleuth, normalize all the studies
# merge into final expression table
# 
