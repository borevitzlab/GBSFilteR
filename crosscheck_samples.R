args <- commandArgs(trailingOnly=T)
setwd("/home/kevin/UniWork/BIOL3157/Borevitz ASC/pel_crosscheck/")
print("filtered")
pel2b.filtered <- read.csv("pel2b.filtered.csv")
pel3.filtered <- read.csv("pel3.filtered.csv")
pel2bn3.filtered <- read.csv("pel2bn3.filtered.csv")
pel2bnpel3.filtered <- c(as.character(pel2b.filtered$tagseqs), as.character(pel3.filtered$tagseqs))
table(!is.na(match(pel2b.filtered$tagseqs, pel3.filtered$tagseqs)))
table(!is.na(match(pel2bn3.filtered$tagseqs, pel2b.filtered$tagseqs)))
table(!is.na(match(pel2bn3.filtered$tagseqs, pel3.filtered$tagseqs)))
table(!is.na(match(pel2bnpel3.filtered, pel2bn3.filtered$tagseqs)))

print("unfiltered")
pel2b.unfiltered <- read.csv("pel2b.fas.csv")
pel3.unfiltered <- read.csv("pel3.fas.csv")
pel2bn3.unfiltered <- read.csv("pel2bn3.fas.csv")
pel2bpel3.unfiltered <- c(as.character(pel2b.unfiltered$seqs), as.character(pel3.unfiltered$seqs))

table(!is.na(match(pel2b.unfiltered$seqs, pel3.unfiltered$seqs)))
table(!is.na(match(pel2bn3.unfiltered$seqs, pel2b.unfiltered$seqs)))
table(!is.na(match(pel2bn3.unfiltered$seqs, pel3.unfiltered$seqs)))
table(!is.na(match(pel2bpel3.unfiltered, pel2bn3.unfiltered$seqs)))

print("key")
key <- read.csv("key.pel2bn3.csv")
pel2b.key <- key[key$PlateName == "Pel2b",]
pel3.key <- key[key$PlateName == "Pel3",]
techrep.samples <- match(as.character(pel2b.key$Sample), as.character(pel3.key$Sample))
print(techrep.samples)


