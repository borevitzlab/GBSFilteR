# Tricia's 7p by Justin on Aug 8,2015
fileLoc <- "Downloads" # "http://gduserv.anu.edu.au/~borevitz/tassel/FOLDERNAME/HapMap/"
fileName <- "Tricia7p.hmc.txt"
setwd(fileLoc)
#hmc <- read.table(paste(fileLoc,fileName,sep="/"),header=T)
hmc <- read.table(fileName,header=T)
#hmc <- read.table("HapMap.hmc.txt",header=T)
# this is taking a while to load, it is big! 332Mb

#exclude standard header columns
std.head <- c("rs","HetCount_allele1","HetCount_allele2","Count_allele1","Count_allele2","Frequency")

#setup matrix of sample rows 1x with allele pairs 2x
hmc.allele <- apply(hmc[,!colnames(hmc)%in%std.head],1,function (x) unlist(strsplit(x,split="|",fixed=T)))
sampN <- nrow(hmc.allele)/2
snpN <- ncol(hmc.allele)
names.list <- list(colnames(hmc[,!colnames(hmc)%in%std.head]),paste(rep(hmc$rs,each=2),1:2,sep="_")  )
#short.names <- matrix(unlist(strsplit(colnames(hmc[,!colnames(hmc)%in%std.head]),split="_")),nr=5)#[1,]

lapply(strsplit(colnames(hmc[,!colnames(hmc)%in%std.head]),split="_"),length)

#plate <- short.names[3,]
#short.names <- paste(short.names[1,],short.names[3,],sep="_")
#names.list <- list(short.names,paste(rep(hmc$rs,each=2),1:2,sep="_")  )

#split genotypes into allele groups 2x samples 1x snps
hmc.allele2 <- matrix(ncol=2*snpN,nrow=sampN,dimnames = names.list)
# fill allele pairs
hmc.allele2[,seq(1,snpN*2,2)] <- as.numeric(hmc.allele[seq(1,sampN*2,2),])
hmc.allele2[,seq(2,snpN*2,2)] <- as.numeric(hmc.allele[seq(2,sampN*2,2),])
save(hmc.allele2,file="Tricia7p.hmc.allele2.RData",compress=T)

#call Presence/Absense
hmc.allele01 <- hmc.allele2
hmc.allele01[hmc.allele2 != 0] <- 1

# look at coverage across samples
reads.samp <- rowSums(hmc.allele2)
alleles.samp <- rowSums(hmc.allele01)
# VISUALLY INSPECT BAD SAMPLES
pdf("Tricia7pReadsAllelesSamp.pdf")
plot(alleles.samp,reads.samp,xlab = "Alleles Called per Sample", ylab = "Total Reads per Sample",main = "Reads vs Alleles per Sample")
s.cuts <- 1500 #need to edit this for each experiment
abline(v=s.cuts)
abline(a=0,b=5,col="red")
abline(a=0,b=10,col="green")
abline(a=0,b=20,col="blue")
dev.off()
keep <- alleles.samp>s.cuts
table(keep)
hmc01 <- hmc.allele01[keep,]

# look at read coverage
# CONSIDER PLATE EFFECT ON ALLELE COUNTS
reads.allele <- colSums(hmc.allele2[keep,])
samps.allele <- colSums(hmc01)
plot(samps.allele, reads.allele,xlab = "Samples with Allele Called", ylab = "Total Reads per Allele",main = "Reads vs Alleles per Locus",
     pch=".",ylim=c(0,10000))
abline(a=0,b=3,col="red")
abline(a=0,b=10,col="green")
abline(a=0,b=40,col="blue")

quantile(reads.allele,1:20/20)
quantile(samps.allele,1:20/20)
freq.allele <- samps.allele/sum(keep)
# remove rare and repeat alleles
s.cuts <- c(-0.1,0.05,0.95,1.1) # must be observed in 5% of samples but not more that 95%
s.cuts <- c(-0.1,0.1,0.9,1.1) # must be observed in 10% of samples but not more that 90%

abline(v=s.cuts*sum(keep))
keep.allele <- cut(freq.allele,s.cuts,labels = c("rare","mid","repeat"))
table(keep.allele)
table(is.na(keep.allele) ) #verify nothing missing because exactly on threshold
hmc01 <- hmc01[,keep.allele=="mid"]  # strange bug, check diminsions
dim(hmc01)

#finally relatedness
# understand 'binary' distance which disregards shared absence
test <- rbind(c(0,0,1,1,1,1,1,1,1,1),
              c(0,0,0,1,1,1,1,1,1,1),
              c(1,1,1,0,0,0,0,0,0,0),
              c(1,1,0,0,0,0,0,0,0,0))
as.matrix(dist(test,method="manhattan"))
as.matrix(dist(test,method="binary"))
plot(hclust(dist(test,method="binary")))
#and read data

dist01 <- dist(hmc01,method="binary")

hc <- hclust(dist01)
pdf(file="Tricia7pDendrogramOdd.pdf",paper="a4r")
plot(hc,cex=0.1)
dev.off()
abline(h=.9)

# cut out odd species as listed below.
cat.tree <- cutree(hc,h=.95)
table(cat.tree)
odd.species <- which(cat.tree > 1)
hmc01 <- hmc01[-odd.species,]
dim(hmc01)

dist01 <- dist(hmc01,method="binary")
hc <- hclust(dist01)
pdf(file="Tricia7pDendrogram.pdf",paper="a4r")
plot(hc,cex=0.1)
dev.off()



which(cat.tree==3)
which(cat.tree==2)
(if i==0){
> which(cat.tree==3)
LBM_1_04_S45649_1_none_A01_X5  LBM_1_08_S45653_1_none_A01_X5 
291                            292 
LBM_1_12_S45657_1_none_A01_X5  LBM_2_03_S45670_1_none_A01_X5 
293                            298 
LBM_2_05_S45672_1_none_A01_X5  LBM_2_07_S45674_1_none_A01_X5 
299                            301 
LBM_5_11_S45679_1_none_A01_X5 LBM_5_13R_S45682_1_none_A01_X5 
302                            304 
LBM_5_15_S45684_1_none_A01_X5  LBM_3_07_S45695_1_none_A01_X5 
306                            309 
M.3_S45051_1_none_A01_X5 
325 
> which(cat.tree==2)
KMCH_01_01_S45528_1_none_A01_X5 KMCH_01_02_S45536_1_none_A01_X5 
69                              70 
KMCH_01_03_S45544_1_none_A01_X5 KMCH_01_04_S45552_1_none_A01_X5 
71                              72 
KMCH_01_05_S45560_1_none_A01_X5              LBM_6_01_merged_X3 
73                             285 
LBM_5_06_S45460_1_none_A01_X5              LBM_6_02_merged_X3 
286                             288 
LBM_6_03_merged_X3   LBM_1_02_S45647_1_none_A01_X5 
289                             290 
LBM_2_01_S45668_1_none_A01_X5   LBM_2_06_S45673_1_none_A01_X5 
296                             300 
LBM_5_12_S45680_1_none_A01_X5   LBM_5_14_S45683_1_none_A01_X5 
303                             305 
LBM_5_17_S45686_1_none_A01_X5    LBM_3_1_S45689_1_none_A01_X5 
307                             308 
LBM_1_25_S45726_1_none_A01_X5   LBM_6_04_S45752_1_none_A01_X5 
310                             311 
LBM_6_05_S45753_1_none_A01_X5   LBM_6_06_S45754_1_none_A01_X5 
312                             313 
#coloring of trees from Aaron Chuah
}


colors <- as.numeric(factor(metadata$Site[match(ind.names,metadata$Tube.label)]))

pl.out <- plot(hclust(as.dist(dist01)),cex=0.5)

library(RColorBrewer)
category.vals<<-as.character(levels(as.factor(categories)))
color.pal<<-c('black',brewer.pal(length(category.vals),"Paired"))

colLab <- function(n) { #helper function to color dendrogram labels
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol<-color.pal[1]
    for(i in 1:length(category.vals)) {
      family<-matrix(unlist(strsplit(a$label,split="-")),nrow=2)[2]
      if(!is.na(pmatch(category.vals[i],family))) {
        labCol<-color.pal[i+1]
      }
    }
    attr(n, "nodePar") <- c(a$nodePar, lab.col=labCol, pch=NA)
  }
  n
}


hc <- hclust(as.dist(dist01))
hc <- dendrapply(as.dendrogram(hc,hang=0.1),colLab)
par(cex=0.3)
plot(hc)
