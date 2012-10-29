
args <- commandArgs(trailingOnly=T)
# 

tassel.result <- read.csv(args[1])
tassel.result <- tassel.result[,2:length(tassel.result)]

names.list <- unlist(strsplit(gsub("^([[:alpha:]]+)\\.", "\\1_", as.character(names(tassel.result)), "."), "_"))
names <- names.list[seq(1,length(names.list),2)]
names.num <- as.numeric(factor(names))

correlations = numeric(length=length(names.num))

for (iii in seq(1,length(names.num))){
  correlations[iii] = cor(as.numeric(factor(names)),as.numeric(tassel.after[iii,]), use = "complete")
}
