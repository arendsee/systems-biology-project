# Get a map from transcription factor name to AGI (TAIR ID)

require(rvest)
require(magrittr)

d <- read_html('http://www.athamap.de/documentation_matrix_based.php') %>%
     html_node('.tableborder table') %>%
     html_table(header=TRUE)

d$Factor <- gsub('.Target.*', '', d$Factor)
d$AGI <- gsub('At(\\d)g', 'AT\\1G', d$AGI)

d <- d[c('Factor', 'Family', 'AGI')]

names(d) <- c('factor', 'family', 'locus')
d <- subset(d, grepl('^AT\\dG\\d{5}$', d$locus, perl=TRUE))

write.table(d, file='OUTPUT/tf.tab', quote=F, sep="\t", row.names=F)
