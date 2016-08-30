### Shuangshuang Lu, Jared Streich, Justin Borevitz 
### code notes, not clean
### Sept 1st, 2016


#change working directory
#setwd("Downloads")

# read in numberic raw file  # CALinesHapMapFilter.R

###################
# this is the genotype data already filtered
#https://github.com/borevitzlab/GBSFilteR/blob/master/LuBrachyCA.txt
geno <- read.table("https://raw.githubusercontent.com/borevitzlab/GBSFilteR/master/LuBrachyCA.txt",header=T,na.strings=".")

str(geno)

#strip marker info
g <- geno[,-(1:11)] 

# KDM HACKING
g2 = g
g2 = sapply(g, function (x) as.numeric(as.character(x)))
image(g2)
g = g2

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
names(g) <- names.mat[1,]
#image raw data
jpeg(file="bd6.jpg",width = 1024, height = 768)
image(1:nrow(g),1:ncol(g),as.matrix(g),xlab = "SNPs", ylab = "samples", 
       main = paste(c("SNPs","samples"),dim(g)) )
dev.off()
## polyploids look obvious in last plate and end of plate 3 ~110 samples

# identify the number of genotype calls by marker
samplecalls <- apply(g,1, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls,breaks=ncol(g))
samp.thres <- 75
abline(v=samp.thres,lty=2)

#test a couple thresholds
table(samplecalls > samp.thres)
#FALSE  TRUE 
#33376 34970 ## original SNPs called
g.snp <- g[samplecalls > samp.thres,]

# look at counts across samples
snpsPerSample<-apply(g.snp,2,function(x) sum(!is.na(x)));
snpsPerSample<-apply(g,2,function(x) sum(!is.na(x)));
pdf("snpsSampleCALines.pdf")
hist(snpsPerSample,breaks = 40)
#test some sample count thresholds
snpCallThresh <- 1600
abline(v = snpCallThresh,lty= 2)
dev.off()
table(snpsPerSample > snpCallThresh)
names(g.snp)[snpsPerSample < snpCallThresh]

names(g.snp)[snpsPerSample < snpCallThresh]
 #[1] "X10bx2s"   "X2sx5m"    "ALM.1"     "BDTR.13o"  "BRO1.1"    "BUR.7"     "BdTR.8c"   #"CFW.3"     "CFWA.1"    "GUL1.3"    "Kah.3x10b" "Koz.4x1f"  "LAL.4"     "LAL5.2out" #"MEN.2029" 
#[16] "MRY.8"     "OSB.2"    


#FALSE  TRUE 
#    5    91 

names(g.snp)[snpsPerSample < snpCallThresh]
[1] "X2sx5m"   "BDTR.13o" "GUL1.3"   "MRY.8"    "OSB.2"   

g.minor <- g.snp[,snpsPerSample > snpCallThresh]

names.mat <- names.mat[,snpsPerSample > snpCallThresh]

# rare varients
# 1 het
rare.var <- rowSums(g.minor,na.rm=T) == 1
table(rare.var)
g.minor <- g.minor[!rare.var,]
# 1 hom
rare.var <- (rowSums(g.minor==2,na.rm=T) == 1) & (rowSums(g.minor,na.rm=T) == 2)
table(rare.var)
g.minor <- g.minor[!rare.var,]

# paralogs?
image(1:nrow(g.minor),1:ncol(g.minor),g.minor==1)
# paralogs/duplicate markers 
het.count <- rowSums(g.minor==1,na.rm=T)
hist(het.count,breaks=ncol(g.minor))
paraThresh <- 40
abline(v = paraThresh, lty = 2)
table(het.count < paraThresh)
#FALSE  TRUE 
#18644 16134 
g.minor <- g.minor[het.count < paraThresh,]
# how many total calls?
table(is.na(g.minor))
image(as.matrix(g.minor))

# minor allele frequency distribution
hist(rowSums(g.minor,na.rm=T)/2)

test <- t(as.matrix(g))
test[test==0] <- "c"
test[test==1] <- NA
test[test==2] <- "t"

require.package(ape)
library(ape)
hap.genotypes.bin <- as.DNAbin(test)
## "raw" is the proportion of sites that differ between each pair of sequences.
hap.gene.dist <- dist.dna(hap.genotypes.bin, model = "raw" , pairwise.deletion=T)

hist(hap.gene.dist,breaks=100)
image(1:nrow(test),1:nrow(test),as.matrix(hap.gene.dist))
write.csv(as.matrix(hap.gene.dist),file = "bd6geneticDist.csv")

pdf(file = "CaliforniaBrachyKinshipMatrix.pdf", paper = "a4", width = 11, height = 15)
image(1:nrow(test),1:nrow(test),as.matrix(hap.gene.dist))
dev.off()


pdf(file = "dendrogramcalibrachy.pdf",paper = "a4",width = 11, height = 15)
plot(nj(hap.gene.dist),cex=0.3)
dev.off()

pdf(file = "caliBrachy_no_hang.pdf",paper = "a4",width = 11, height = 15)
plot(hclust(hap.gene.dist),cex=0.3)
dev.off()

swil.1 <- cmdscale(hap.gene.dist, k = 8)

swil.1_df <- data.frame(swil.1)
tot.var <- sum(sapply(swil.1_df, var))
per.var <- round((sapply(swil.1_df, var)/tot.var)*100,2)

pdf(file = "caliPCoA.pdf",paper = "a4r",width = 11, height = 15)
plot(swil.1[,2]~swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachy PCoA" )
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels =rownames(swil.1) ,
                            cex=0.4,col = names.mat[3,] )
#legend(500,400,paste("plate",c(2,3,4,6)),col = c(2,3,4,6),lty=1)
dev.off()

pdf(file = "caliBrachyPCoAsm.pdf",paper = "a4r",width = 11, height = 15)
plot(swil.1[,2]~swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachy PCoA",type="n" )
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels =rownames(swil.1) ,
                            cex=0.6,col = names.mat[3,] )
legend(500,400,paste("plate",c(2,3,4,6)),col = c(2,3,4,6),lty=1)
dev.off()


### The other way to calculate distance
snp.cor.dist <- as.dist(1-cor(g.minor,use = "pairwise.complete.obs")^2)

pcoa.out <- pcoa(snp.cor.dist)
biplot(pcoa.out,col=rep(gps$color[index],2),cex=.1)

pdf(file = "calicorhclust.pdf",paper = "a4r",width = 11, height = 15)
plot(hclust(snp.cor.dist),cex = 0.3)
dev.off()

########################## CABrachy.popgen.short2.R  

library("ape") 
library("pegas") 
library("seqinr") 
library("ggplot2")
library("adegenet") 

setwd("/Users/Sabrina/Documents/*R/0.MyRscript/adegenet.Practice.Sept2015")

# read in numberic raw file
geno <- read.table("P9_10_no1.txt",header=T,na.strings=".")

#strip marker info
g <- geno[,-(1:11)] 
dim(g) #dim 5516x176
g1<-g[,9:174]
g1[1:2,1:10]
colnames(g1)<-pop
dim(g1)
write.table(g1,file="CAbrachy_treemix.txt")
rownames(g1) 
S1<-ddply(g1,.(colnames(g1)),summarise,zeros=sum(rownames=="0"))

# switch columns and rows (transpose) 
tg<-data.frame(t(g))
#subset tg, remove non-CA accessions
tg1<-tg[9:174,] #dim 166x5516
#remove rownames of tg1

#tried for treemix
g1<-data.frame(t(tg1))
dim(g1)
colnames(g1)
pop.t<-data.frame(t(pop))
##
colnames(tg1)[0]<-"Pop"
tg2<-tg1
rownames(tg2)<-c()
CA.brachy.SNPdf<-cbind(acc,pop.int,tg2)


write.table(CA.brachy.SNPdf,file="CA.brachy.SNPdf.txt",sep=" ",row.names = FALSE)

CA.brachy.SNPdf<-read.table("CA.brachy.SNPdf.txt",header=T,na.strings="NA")

CA.brachy.SNPdf[is.na(CA.brachy.SNPdf)]<-"-9" #replace NA with -9 for STRUCTURE
write.table(CA.brachy.SNPdf,file="CA.brachy.SNPdf.txt",sep=" ",row.names = FALSE)

snpdf<-read.table("CA.brachy.SNPdf.txt")
snpdf[1:2,1:10]
snpdf$V2
colnames(snpdf)

#count numbers of zeros and twos at each loci for each population
# "paste" was used to merge zero and two into one column for each loci
snpsum<-aggregate(snpdf[paste0('V',3:5518)],by=snpdf['V2'],
                  function(x) paste(zero=sum(x=="C/C"),two=sum(x=="A/A"),sep=','))

snpsum[1:10,1:10]
rownames(snpsum)<-snpsum[,1]
colnames(snpsum)<-NULL
snpsum<-snpsum[,-1]
write.table(snpsum,file="snpsum.txt")

snpsum<-read.table("snpsum.txt")
dim(snpsum)
tsnpsum<-t(snpsum) # transpose data to format for treemix
dim(tsnpsum)
tsnpsum[1:10,1:2]
rownames(tsnpsum)<-NULL
pop.name<-c("Fremont1","Fremont2","KirkCreek","MorroBay1","MorroBay2","SantaBarbara","LosAngeles","LongBeach","DanaPoint","SanDiego","ChulaVista","DonPedroReservior1","DonPedroReservior2","Winters1","Winters2","Folsom1","Folsom2","Folsom3","Oroville","Chico")
colnames(tsnpsum)<-pop.name
write.table(tsnpsum,file="CAbrachy_treemix.txt",row.names=F)
#test :
tree<-read.table("CAbrachy_treemix.txt",header=T)
dim(tree)
tree[1:2,1:10]
#after running treemix from terminal, plotting via R:
source("/Users/Sabrina/Documents/Software/Treemix/treemix-1.12/src/plotting_funcs.R")
setwd("/Users/Sabrina/Documents/Software/Treemix/brachydata/")
pdf(file="CAbrachy-Tree-SantaBarbara.pdf",width=12,height=10)
plot_tree("CAbrachy")
dev.off()

# Generate "genlight" object from SNP data  
x<-new("genlight",tg1)
str(x)
write.table(x,file="genlight_tg1.txt")

l=read.table("genlight_tg1.txt") # try to save genlight data and read it from folder
# simplify accession names and add population info
acc<-read.table("indi.names.loc_CA.txt")
pop<-read.table("pop_CA.txt")
pop.int<-read.table("pop_integer.txt")
acc<-as.character(acc[,1])
pop<-as.character(pop[,1])

indNames(x)<-acc #add accession (individual) name
pop(x)<-pop # Assign Pop information
#ploidy(x)<-2 # try set ploidy level to 2
x2=x
ploidy(x2)=2
# principal component analysis ####
pca1<-glPca(x) # Ask to select the number of axes. Select 2?
scatter(pca1,posi="bottomright") #

#neighbour-joining tree:
library(ape)
tre<-nj(dist(as.matrix(x)))
plot(tre,typ="fan",cex=0.7)
title("NJ tree of CA Brachy")

# use colorplot to combine PCA and NJ graphs####
pdf("CAbrachy_PCcolorplot.pdf",width=10,height=8)
myCol<-colorplot(pca1$scores,pca1$scores,transp=TRUE,cex=4)
#text(pca1$scores[,1],pca1$scores[,2],x$pop,cex=1)
abline(h=0,v=0,col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2,posi="topright",inset=.05,ratio=.3)
dev.off()

pdf("CAbrachy_NJtree.pdf",width=20,height=20)
plot(tre,typ="fan",show.tip=TRUE)
tiplabels(pch=20,col=myCol,cex=4)
title("NJ tree of CA Brachypodium")
dev.off()

# Discriminant Analysis of Principal Components (DAPC) ####
# first use find.clusters to transform genlight object via PCA
grp<-find.clusters(x,max.n.clust=40) 
# number of PCs retained: 200 (retain all information here)
#number of clusters chosen: 9
names(grp)
head(grp$Kstat,8)
head(grp$grp,10)
grp$size
grp$grp

# the DAPC
dapc11 <- dapc(H3N2, pop=H3N2$other$epid, n.pca=30,n.da=6)
dapc1<-dapc(x,pop=x$other$epid,n.pca=20)
dapc1<-dapc(x,grp$grp) #retained 20 PCs
#the result displays a barplot of eigenvalues for the discriminant analysis
# 8 discriminant functions were saved
scatter(dapc1, posi.da="topleft", bg="white", pch=17:22) 

#scatter plot with 9 clusters####
pdf("CAbrachy_9clusters.pdf",width=20,height=20)
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=3,  mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",posi.leg = "bottomleft", leg=TRUE, txt.leg=paste("Cluster",1:9))
text(pca1$scores[,1],pca1$scores[,2],x$pop,cex=1)
# if want to remove label of clusters: use clab=0
par(xpd=TRUE) 
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=1, lwd=5, col="black")

myInset<-function(){
  temp<-dapc1$pca.eig
  temp<-100*cumsum(temp)/sum(temp)
  plot(temp,col=rep(c("black","lightgrey"),c(dapc1$n.pca,1000)),ylim=c(0,100),
       xlab="PCA axis",ylab="Cumulated variance (%)",
       cex=1,pch=20,type="h",lwd=2)
}
add.scatter(myInset(),posi="bottomright",inset=c(-0.05,-0.01),ratio=.2,bg=transp("white"))
dev.off()

#plotting the densities of individuals on a given discriminant function with different colors for different groups
#discriminant function 1
scatter(dapc1,1,1, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

set.seed(4)
# show the most contributing alleles: contributions above a given threshod (.005)
contrib<-loadingplot(dapc1$var.contr,axis=2,thres=.005,lab.jitter=1)


#stAMPP for F-stats and genomic matrix ####
library(StAMPP)
# use genlight object "x" 
x1<-x
ploidy(x)<-1
x1
x2.freq<-stamppConvert(x2,type="genlight") 
x.freq<-stamppConvert(x,type="genlight") # not sure whether the conversion of 0,1,2 SNP data into 0.25, 0.5, 0.75 is correct?

# calculate Nei's Genetic Distance by populations ####
Nei.pop<-stamppNeisD(x.freq,pop=TRUE)
Nei.pop2<-stamppNeisD(x2.freq,pop=TRUE) #ploidy=2
write.csv(Nei.pop,"CAbrachy_pop_Nei_hap.csv") #ploidy=1
# calculate Nei's Genetic Distance by individuals
Nei.ind<-stamppNeisD(x.freq,pop=FALSE)
#write.csv(Nei.ind,"CAbrachy_ind_Nei_distance.csv")
#Amova:
amo=stamppAmova(Nei.ind,x.freq,perm=100)
summary(amo)
amo$varcomp

# Calculate Fst####
fst<-stamppFst(x.freq,nboots=100,percent=95,nclusters=2) # probably should increase nboots #diploid here
fst2<-stamppFst(x2.freq,nboots=100,percent=95,nclusters=2) 
fst3<-stamppFst(x2.freq,nboots=1000,percent=95,nclusters=2) #increase nboots
fst$Fsts
fst3$Fsts
identical(fst$Fsts,fst1$Fsts) # Results are the same for treating as haploids or diploids
fst3$Pvalues
capture.output(fst$Fsts,file="CABrachy.fst.csv")

capture.output(fst2$Fsts,file="CABrachy.fst.csv")
capture.output(fst3$Fsts,file="CABrachy.1000bootfst.csv")
write.csv(fst3$Fsts,file="CABrachy.1000bootfst.csv")

capture.output(fst2$Pvalues,file="CABrachy.fstPvalues.csv")
# Calculating a genomic relationship matrix ####
x.g<-stamppGmatrix(x.freq)

# Generate a heatmap
pdf("CAbrachy_heatmap_no1.loc-trash.pdf",width=25,height=25)
heatmap(x.g,
        margins=c(5,10),
        xlab="Assigned Individual ID (bottom right of x-axis labels corresponds to top right of y-axis labels ",ylab="Individual Name",main="Genetic Clusters of CA Brachypodium Populations (cells from yellow to red colors indicate increasing genetic distance)")
dev.off()

library(gplots)
colors <- rev(heat.colors(256))
colbr <- c(seq(-1, 1, len=length(colors)+1))
heatmap.2(x.g,trace='none',
          margins=c(5,10),
          xlab="Assigned Individual ID (bottom right of x-axis labels corresponds to top right of y-axis labels ",ylab="Individual Name",main="Genetic Clusters of CA Brachypodium Populations (cells from yellow to red colors indicate increasing genetic distance)")

