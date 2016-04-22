# # Install rhdf5 with the following code, if needed
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")

# # Install sleuth with the following if needed
# source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")    # only if devtools not yet installed
# biocLite("pachterlab/sleuth")

require("sleuth")
require('magrittr')
require('rhdf5')

base_dir=file.path(getwd(), 'output')

d <- read.delim('data/conditions.tab')
d <- d[,c('run_accession', 'condition', 'drop')] 
d <- subset(d, drop == 0)[c(1,2)] %>% droplevels
names(d) <- c('sample', 'condition')
d$path <- sapply(d$sample, function(id) file.path(base_dir, id))

# For some reason, this entry, the drought condition, died while running kallisto
d <- subset(d, sample != 'SRR932112') %>% droplevels

so <- sleuth_prep(d, ~ condition)

# Memory death!
so <- sleuth_fit(so)
