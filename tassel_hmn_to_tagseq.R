args <- commandArgs(trailingOnly=T)
# args[1] = "pel2b_crosscheck.hmn.txt"
# args[2] = "pel2b_crosscheck.fas.csv"
# args[3] = "pel2b_crosscheck"

### READ DATA
tassel <- read.table(args[1], header=T)
tassel.seqs <- read.csv(args[2])

###MATCH ROWS
seq.tp.match <- match(as.character(tassel.seqs$seq.names), tassel[,1])
tassel[,1] <- tassel.seqs[seq.tp.match,3]

#write result
write.csv(tassel, file=paste(args[3], "hmnfas.csv", sep="."))