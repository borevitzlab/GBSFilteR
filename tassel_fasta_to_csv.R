library(seqinr)

args <- commandArgs(trailingOnly=T)
fasta <- read.fasta("HapMap.fas.txt", as.string=T, seqonly=F, forceDNAtolower=F)

seq.names.raw <- as.character(matrix(unlist(names(fasta))))
seq.names.raw <- seq.names.raw[seq(1,length(seq.names.raw),2)]
seq.names <- as.character(matrix(unlist(strsplit(seq.names.raw, "_"))))
seq.names <- seq.names[seq(1,length(seq.names),2)]

seqs <- as.character(matrix(unlist(fasta)))
seqs <- seqs[seq(1,length(seqs),2)]

table <- data.frame(seq.names=seq.names, seqs=seqs)

write.csv(table, file="HapMap.fas.csv")

