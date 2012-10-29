args <- commandArgs(trailingOnly=T)
#1: Tassel hmn.txt table
#2: fasta csv
#3: hmnfas tagseq output file

args[1] = "pel2b_crosscheck.hmn.txt"
args[2] = "pel2b_crosscheck.fas.csv"
args[3] = "pel2b_crosscheck.hmnfas.csv"

tassel <- read.table(args[1], header=T)
tassel.seqs <- read.csv(args[2])
seq.tp.match <- match(as.character(tassel.seqs$seq.names), tassel[,1])
tassel[,1] <- tassel.seqs[seq.tp.match,3]

write.csv(tassel, file=args[3])