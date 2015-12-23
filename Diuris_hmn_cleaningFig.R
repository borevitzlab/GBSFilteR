# Diuris .hmn file cleaning of data
geno <- read.table(file.choose(), head = TRUE, na.strings=".",row.names = 1)

geno <- geno[,-(1:10)] 

#sort columns
geno <- geno[,order(names(geno))]


#remove alpine individuals
g.m <- g.minor[,-(1:15)] 

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(geno),split="_")),nrow=6)
names(geno) <- names.mat[1,]

#image raw data

jpeg(file=paste0(dataset,"_snps.jpg"),width = 1024, height = 768)

image(1:nrow(geno),1:ncol(geno),as.matrix(geno),xlab = "SNPs", ylab = "samples", main = paste(c("SNPs","samples"),dim(geno)) )

dev.off()


# look at counts across samples

snpsPerSample<-apply(geno,2,function(x) sum(!is.na(x)));

hist(snpsPerSample,breaks = 40)


#test some sample count thresholds

snpCallThresh <- 4000

abline(v = snpCallThresh,lty= 2,col="orange")

table(snpsPerSample > snpCallThresh)

g.minor <- geno[,snpsPerSample > snpCallThresh]

names.mat <- names.mat[,snpsPerSample > snpCallThresh]

names(g.minor) <- names.mat[1,]


# identify the number of genotype calls by marker

samplecalls<-apply(g.minor,1, function(x) sum(!is.na(x)));


#look at distribution and try to determine a threshold for markers

jpeg(file=paste0(dataset,"histCounts.jpg"),width = 1024, height = 768)

hist(samplecalls,breaks=ncol(geno))

# stringent data set = 200; relaxed data set = 70
samp.thres <- 200

abline(v=samp.thres,lty=2)

dev.off()


#test a couple thresholds

table(samplecalls > samp.thres)

g.snp <- g.minor[samplecalls > samp.thres,]




# rare variants

# 1 het

rare.var <- rowSums(g.snp,na.rm=T) == 1

table(rare.var)

g.minor <- g.snp[!rare.var,]

# 1 hom

rare.var <- (rowSums(g.snp==2,na.rm=T) == 1) & (rowSums(g.minor,na.rm=T) == 2)

table(rare.var)

g.minor <- g.snp[!rare.var,]

dim(g.minor)



# paralogs?

image(g.snp==1)

# paralogs/duplicate markers

het.count <- rowSums(g.snp==1,na.rm=T)

hist(het.count)

paraThresh <- 45

abline(v = paraThresh, lty = 2)

table(het.count < paraThresh)


g.snp <- g.snp[het.count < paraThresh,]

# how many total calls?

table(is.na(g.snp))

image(as.matrix(g.snp))


# minor allele frequency distribution

hist(rowSums(g.snp,na.rm=T)/2)


# distance matrices

snp.dist <- dist(t(g.snp))

write.table(g.snp, file="fil_200_hmn.txt")
library("ape")

pcoa.out <- pcoa(snp.dist)

pdf(paste0(dataset,"PCOAdendrogram.pdf"),paper = "a4r",width = 11, height = 15)
par(cex=0.8)
biplot(pcoa.out)
plot(hclust(snp.dist),cex=0.8)

# Print pcoa scores
ax1<-pcoa.out$vectors[,1]
melt(ax1)
#axis 2
ax2<-pcoa.out$vectors[,2]
melt(ax2)

################# Fig4 
library(adegenet)
m<-read.table(file.choose(), row.names=1)

# convert to genind obj
N<-df2genind(m, missing='NA', sep="|")


#read in population table
t<-read.table(file.choose(), header=TRUE, row.names=1)
pop(N)<-t[,1]


col17=c('black','yellow','blue','blue','blue','blue','red','blue','green','green','green','yellow','green','green','red','red','blue')


library(ape)
#PCOA
X <- scaleGen(N, scale=FALSE, miss="mean") 
p <- dudi.pco(dist(X), scannf=FALSE, nf=3)
plot(as.phylo(hclust(dist(N))), type="unrooted", tip.color=col17[pop(N)])




#Building tree from Fst Values
fst<-pairwise.fst(N,res.type="matrix")

#delete pops with missing data
fst2<-fst[-12,-12]
fst2<-fst2[-11,-11]
fst2<-fst2[-8,-8]
fst2<-fst2[-5,-5]
fst2<-fst2[-3,-3]
fst2<-fst2[-1,-1]

bin.tree <- nj(fst2) 
plot(bin.tree, type="unr", font=2) 
add.scale.bar()
