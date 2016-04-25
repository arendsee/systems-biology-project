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

stuid1 = c(
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

stuid3 = c(
  # Developmental stages
  'ERP004034', # inflorecence meristem, floral meristem and flower 
  'ERP005391', # PAIRED: three stages of somatic embryo development
  'SRP007113', # female gametophyte
  'SRP009369', # PAIRED: root hair 
  'SRP035269', # seed development (3 lines) X (3 dev stages) X (3 bioreps)
  # hormone perturbations
  'SRP010642', # brassinosteroid and gibberillin treatments, 
  # abiotic stresses
  'SRP018404', # root response to nitrate
  'SRP003864', # copper defficiency (WT and spl7-KO) X (normal, low [Cu])
  'SRP027256', # drought condition versus WT
  'SRP044814', # Fe stress
  # biotic stress
  'SRP034715', # Response to mildew (huge)
  'SRP063017'  # Nematode infection
)


fjoin <- function(x) sprintf("'%s'", paste0(x, collapse="', '"))

getRuns <- function(s){
  sprintf("'%s'", paste0(s, collapse="', '")) %>%
    sprintf(fmt='select * from sra_ft where study_accession in ( %s )') %>%
    query()
}

runs1 <- getRuns(stuid1)
runs2 <- getRuns(stuid2)
runs3 <- getRuns(stuid3)

out <- runs3[c(
  'run_accession',
  'experiment_accession',
  'sample_accession',
  'submission_accession',
  'study_accession',
  'sample_name',
  'sample_alias',
  'experiment_name',
  'experiment_alias',
  'experiment_title',
  'library_layout',
  'platform',
  'library_strategy',
  'library_selection',
  'study_name',
  'study_title',
  'study_abstract',
  'study_type',
  'study_description'
)]
write.table(out, file='z.tab', quote=F, sep="\t", row.names=F)
write.table(runs3, file='zfull.tab', quote=F, sep="\t", row.names=F)

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
