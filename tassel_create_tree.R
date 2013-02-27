library(phangorn)

#change working directory
args <- commandArgs(trailingOnly=T)
# args[1] = "pel2b.filtered.csv"
# args[2] = "pel2b"

### READ DATA
# Read genotype data created tassel_filter_hapmap.R
genotype.seqs <- read.csv(args[1])
# remove tag seqs for tree generation
genotype <- genotype.seqs[,seq(2,length(genotype.seqs)-1)] 
row.names(genotype) <- genotype.seqs$X # restore row names


### MAKE TREE
# Make tree, define colour groups
tree <- hclust(as.dist(1-cor(genotype,use = "pairwise.complete.obs")))
tree.groups <- cutree(tree, k=5)
write.csv(tree.groups, file=paste(args[2], "phylogroups.csv", sep="."))

# Assigns each phylogenetic group a colour
plot.col <- rainbow(max(as.numeric(tree.groups)))[as.numeric(tree.groups)]


### MAKE TREE FIGURE
pdf(file=paste(args[2], "tree.pdf", sep="."))
# Plot tree, including coloured boxes around different phylogenetic groups
plot(
  tree, 
  cex = 0.7, #cex scales text by 0.5
  #sub=paste("sample_cutoff=", sample.cutoff, "snp_cutoff=", snp.cutoff),
  xlab=""
)
rect.hclust(
  tree, 
  k=5, 
  border=rainbow(max(as.numeric(tree.groups))), 
  cluster=tree.groups
)

### MAKE BOOTSTRAPPED TREE

#make distance matrix and inital tree with upgma

#snp.cor.dist <- as.dist(1 - cor(genotype,use = "pairwise.complete.obs"))
snp.euc.dist <- dist(t(genotype)) # problems with missing data
initial.tree <- upgma(snp.euc.dist)

## Calculate boostrapped tree
bs <- list()
for (i in 1:100){
  # resample data
  bs.data <- sample(nrow(genotype),replace=T)
  # get dist matrix of resampled data
#  snp.cor.bs <- as.dist(1-cor(genotype[bs.data,],use = "pairwise.complete.obs"))
  snp.cor.bs <- dist(t(genotype[bs.data,],use = "pairwise.complete.obs"))
  # store tree of resampled data
  bs[[i]] <- upgma(snp.cor.bs)
  
  print(i)
}

## Make bootstrapped tree
bs.tree <- plotBS(initial.tree, bs)

## Print Bootstrapped tree
pdf(file=paste(args[2], "treeBS.pdf", sep="."))
plot(initial.tree, cex=0.5)
nodelabels(bs.tree$node.label,frame="none", cex=0.5)
dev.off()
