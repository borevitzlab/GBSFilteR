args <- commandArgs(trailingOnly=T)

pel2b <- read.csv("pel2b.filtered.csv")
pel3 <- read.csv("pel3.filtered.csv")
pel2bn3 <- read.csv("pel2bn3.filtered.csv")

techrep.rows <- c()
match(names(pel2b), names(pel3))
table(!is.na(match(pel2b$tagseqs, pel3$tagseqs)))
table(!is.na(match(pel2bn3$tagseqs, pel2b$tagseqs)))
table(!is.na(match(pel2bn3$tagseqs, pel3$tagseqs)))