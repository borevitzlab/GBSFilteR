#change working directory
args <- commandArgs(trailingOnly=T)
#1: hmnfas file (hmn file with tag sequences)
#2: output file prefix

args[1] = "pel2b_crosscheck.filtered.csv"
args[2] = "pel2b_crosscheck"


# Read genotype data created tassel_filter_hapmap.R
genotype.seqs <- read.csv(args[1])
genotype <- genotype.seqs[,seq(2,length(genotype.seqs)-1)] # remove seqs for tree generation
row.names(genotype) <- genotype.seqs$X # restore row names


### Make tree, define colour groups
tree <- hclust(as.dist(1-cor(genotype,use = "pairwise.complete.obs")))
tree.groups <- cutree(tree, k=5)

# Assigns each phylogenetic group a colour
plot.col <- rainbow(max(as.numeric(tree.groups)))[as.numeric(tree.groups)]


### Export data to gps file
# Creates a data.frame of each sample
# Points are coloured based on cutree phylogenetic groups
# Only samples which passed filtering are used
gps.file <- data.frame(
  Names=names(g.output),
  color=plot.col,
  Lat=names.matrix[8,],
  Long=names.matrix[9,]
)

# Write gps data frame to csv
write.csv(gps.file,file = paste(args[3],"gps.csv", sep="."))
