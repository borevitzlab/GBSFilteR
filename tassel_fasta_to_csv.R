library(seqinr)

args <- commandArgs(trailingOnly=T)
args[1] <- "pel2b_crosscheck.fas.txt"
args[2] <- "pel2b_crosscheck"

# Read data
fasta <- read.fasta(args[1], as.string=T, seqonly=F, forceDNAtolower=F)

seq.names.raw <- as.character(matrix(unlist(names(fasta))))
seq.names.raw <- seq.names.raw[seq(1,length(seq.names.raw),2)]
seq.names <- as.character(matrix(unlist(strsplit(seq.names.raw, "_"))))
seq.names <- seq.names[seq(1,length(seq.names),2)]

seqs <- as.character(matrix(unlist(fasta)))
seqs <- seqs[seq(1,length(seqs),2)]

table <- data.frame(seq.names=seq.names, seqs=seqs)

write.csv(table, file=paste(args[2], "fas.csv", sep="."))

