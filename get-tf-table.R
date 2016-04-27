# Get a map from transcription factor name to AGI (TAIR ID)

require(rvest)
require(magrittr)

d <- read_html('http://www.athamap.de/documentation_matrix_based.php') %>%
     html_node('.tableborder table') %>%
     html_table(header=TRUE)

d$Factor <- gsub('.Target.*', '', d$Factor)
d$AGI <- gsub('At(\\d)g', 'AT\\1G', d$AGI)

d <- d[c('Factor', 'AGI')]

names(d) <- c('factor', 'locus')
d <- subset(d, locus != '-')

write.table(d, file='OUTPUT/factor2agi.tab', quote=F, sep="\t", row.names=F)
