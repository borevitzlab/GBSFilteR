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

### EXPORT FIGURES
pdf(file=paste(args[2], "tree.pdf", sep="."))
# Plot tree, including coloured boxes around different phylogenetic groups
plot(
  tree, 
  cex = 0.7, #cex scales text by 0.5
  sub=paste("sample_cutoff=", sample.cutoff, "snp_cutoff=", snp.cutoff),
  xlab=""
)  
rect.hclust(
  tree, 
  k=5, 
  border=rainbow(max(as.numeric(tree.groups))), 
  cluster=tree.groups
)
dev.off()
