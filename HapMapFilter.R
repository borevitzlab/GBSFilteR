#change working directory
#setwd("Downloads/hapmap2n3")

# read in numberic raw file
geno <- read.table("HapMap.hmn.txt",header=T,na.strings=".")

#strip marker info
g <- geno[,-(1:11)] 
#image raw data
jpeg(file="raw.jpg")
image(as.matrix(g),main = dim(g))
dev.off()

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
# well and plate id in names.mat
names(g) <- names.mat[1,]

# identify the number of genotypes called per sample
sites <- apply(g,1, function(x) sum(!is.na(x)));
hist(sites,breaks=96)
snp.threshold <- 30 ## this depends on sample layout and quality
hist(sites[sites > snp.threshold],breaks= (96-snp.threshold))
g.snp <- g[sites > snp.threshold,]

# drop bad samples with few markers
samples<-apply(g.snp,2,function(x) sum(!is.na(x)))
hist(samples,breaks = 100)

samp.thresh <- 400 # this matters depending on the diversity among samples
table(samples < samp.thresh)
gg <- g.snp[,samples > samp.thresh] #select top 80 samples
# how many total calls?
table(!is.na(g))
#  FALSE    TRUE 
#5195135  474817 
table(!is.na(gg))
# FALSE   TRUE 
#139363 292241 
image(as.matrix(gg))

# paralogs?
image(gg==1)
snp.het.count <- rowSums(gg==1,na.rm=T)
#frequency distribution across samples
hist(snp.het.count, breaks = ncol(gg))
#threshold
paralog.thresh <- 100 # depends on repetitivness 
table(snp.het.count > paralog.thresh)
g.minor <- gg[snp.het.count > paralog.thres,]
write.csv(g.minor,file = "geno-paralog.csv")

snp.cor.dist <- as.dist(1-cor(g.minor,use = "pairwise.complete.obs"))
tree <- upgma(snp.cor.dist)
pdf("tree.pdf")
plot(tree)
bs <- list()
for (i in 1:100){
            bs.mark <- sample(nrow(g.minor),replace=T)
	snp.cor.bs <- as.dist(1-cor(g.minor[bs.mark,],use = "pairwise.complete.obs"))
bs[[i]] <- upgma(snp.cor.bs)
cat(i)
}
install.packages("phangorn",repos = "http://cran.csiro.au/")
library(phangorn)
bstree <- plotBS(tree, bs)
dev.off()
write.nexus(bstree, file = "boot.nex")

