

### EXPORT FIGURES
pdf(file=paste(args[2], "results.pdf", sep="."))
# Minor Allele Frequency histogram of final data
hist(final.min.allele.freq,breaks=50, main="Minor Allele Frequencies")
# Genotype call image of final data
image(
  as.matrix(g.final), 
  main="Genotype Call Data", 
  sub=paste("sample_cutoff=", sample.cutoff, "snp_cutoff=", snp.cutoff),
  col=rainbow(3), 
  ylab="Samples", 
  xlab="SNPs"
)
# Plot tree, including coloured boxes around different phylogenetic groups
plot(
  tree, 
  cex = 0.5, #cex scales text by 0.5
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
