# Tricia's 7p by Justin on Sep 10,2015
fileLoc <- "Downloads" # "http://gduserv.anu.edu.au/~borevitz/tassel/FOLDERNAME/HapMap/"
fileName <- "rerunfinalAug2015.hmc.txt"
fileName <- "http://gduserv.anu.edu.au/~msupple/tricia/final/allmat.HapMap.hmc.txt"
#fileName <- "http://gduserv.anu.edu.au/~msupple/tricia/final/ybmatseed.HapMap.hmc.txt"
#setwd(fileLoc)
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

short.names <- matrix(unlist(strsplit(colnames(hmc[,!colnames(hmc)%in%std.head]),split="_")),nr=7)
names.list <- list(paste(short.names[1,],short.names[2,],sep="-"),paste(rep(hmc$rs,each=2),1:2,sep="_")  )

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
s.cuts <- 2000 #need to edit this for each experiment
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

## need to filter alleles that are present in all lanes!!!
lanes <- metadata$lane[idZ][keep]
lane.count <- matrix(ncol = ncol(hmc01),nr=length(table(lanes)))
for (i in 1:ncol(hmc01)){
lane.count[,i] <- tapply(hmc01[,i],lanes,sum)
} # takes a while 1min

allele.lane <- colSums(matrix(as.logical(lane.count),nr=3))
table(allele.lane)

plot(density(samps.allele ))
lines(density(samps.allele [common.allele]),col="green")
common.allele <- allele.lane == 3
# remove rare and repeat alleles
s.cuts <- c(-0.1,0.05,0.95,1.1) # must be observed in 5% of samples but not more that 95%
s.cuts <- c(-0.1,0.1,0.9,1.1) # must be observed in 10% of samples but not more that 90%

abline(v=s.cuts*sum(keep))
keep.allele <- cut(freq.allele,s.cuts,labels = c("rare","mid","repeat"))
table(keep.allele)
table(is.na(keep.allele) ) #verify nothing missing because exactly on threshold
table(keep.allele,common.allele)

#common.allele
#keep.allele  FALSE   TRUE
#rare   153486  11720
#mid     10301   7428
#repeat      0    273
hmc01 <- hmc01[,keep.allele=="mid" & common.allele]  # strange bug, check diminsions
#hmc01 <- hmc01[,common.allele]  # strange bug, check diminsions
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

# more examples of multiple regression of distance matrices
mat <- rep(LETTERS[1:10], each = 3)
pop <- rep(letters[1:3],each = 10)
kin.mat <- outer(mat,mat, FUN = "==")
pop.mat <- outer(pop,pop, FUN = "==")
# noise + family and pop
gen.mat <- rnorm(30^2) + kin.mat + pop.mat

library(ecodist)
mantel(as.dist(gen.mat) ~ as.dist(kin.mat) + as.dist(pop.mat))
mantel(as.dist(gen.mat) ~ as.dist(pop.mat) )
mantel(as.dist(gen.mat) ~ as.dist(kin.mat) )

MRM(as.dist(gen.mat) ~ as.dist(kin.mat) + as.dist(pop.mat))
MRM(as.dist(gen.mat) ~ as.dist(pop.mat) )
MRM(as.dist(gen.mat) ~ as.dist(kin.mat) )

# principal components first runs the distance matrix, we will do seperately
prcomp.euc <- prcomp(hmc01)

# back to the real data
dist01 <- dist(hmc01,method="binary")
hc <- hclust(dist01)
hist(dist01,breaks=100) # notics 2 peaks for within and between family

image(as.matrix(dist01)[hc$ord,hc$ord])


pc <- cmdscale(dist01,k=30)
var.comp <- apply(pc,2,var)
var.comp <- round(var.comp/sum(var.comp)*100)
colors <- as.numeric(factor(as.character(metadata$PlateName[idx])))
colors <- as.numeric(factor(as.character(metadata$lane[idx])))

plot(pc,xlab = paste("pc1",var.comp[1],"%"),ylab=paste("pc2",var.comp[2],"%"),col = colors)

sampname[order(pc[,1])]
pop[order(pc[,1])]

pdf(file="Tricia683Dendrogram7k.pdf",width=420/25,height=297/25)
plot(hc,cex=0.2,labels = sampname,color = levels(pop))
dev.off()

image(1:length(hc$ord),1:length(hc$ord),as.matrix(dist01)[hc$order,hc$order],zlim=c(0.6,1))
axis(side=1,labels=sampname[hc$order],padj=T,at = 1:length(sampname),cex=0.1)

image(1:length(hc$ord),1:length(hc$ord),as.matrix(dist01)[order(sampname),order(sampname)],zlim=c(0.6,1))
mtext(sampname,side=1, padj = 1,cex=0.1)
dev.off()
abline(h=.9)

temp <- colSums(as.matrix(dist01))
plot(reads.samp[keep],temp)
# cut out odd species as listed below.
cat.tree <- cutree(hc,h=.95)
table(cat.tree)
odd.species <- which(cat.tree > 1)
hmc01 <- hmc01[-odd.species,]
dim(hmc01)

dist01 <- dist(hmc01,method="binary")
hc <- hclust(dist01)
pdf(file="Tricia7pDendrogram.pdf",)
plot(hc,cex=0.1)
dev.off()

metadata <- read.csv("https://github.com/LaMariposa/eucalyptus_data/blob/master/Emelliodora_PlantsSamples.csv")
metadata <- read.csv("Emelliodora_PlantsSamples.csv")

metadata$lane <- "lane"
metadata$lane[metadata$PlateName %in% c("Emel003A","Emel003B")] <- "lane1"
metadata$lane[metadata$PlateName %in% ""] <- "laneX"
metadata$lane[metadata$PlateName %in% paste("Emelliodora_LB_Block",c("09","10","11","12","13"),sep="")] <- "lane2"
metadata$lane[metadata$PlateName %in% paste("Emelliodora_Plate",c(5,6,9,10),"_2015",sep="")] <- "lane3"
metadata$lane[metadata$PlateName %in% paste("Emelliodora-LB-Block",1:4,sep="")] <- "lane4"


#metadata <- read.csv("meta_MS_FP_alt_GPS.csv")

metadata2 <- read.csv("melliodora_lbpops.csv")
idy <- match(as.character(pop),metadata2$pop.name)
long <- metadata2$Long[idy]
lat <- metadata2$Lat[idy]

library(fossil)
geo.dist <- earth.dist(cbind(long,lat))
#make sure all sampleid's have metadata
table(rownames(hmc01)%in%metadata$SampleID)

idx <- match(rownames(hmc01),metadata$SampleID)
idZ <- match(short.names[1,],metadata$SampleID)

sort(as.character(metadata$SampleName[idZ]))

temp <- strsplit(as.character(metadata$SampleName),split = "-")
metadata$family <- unlist(lapply(temp,function (x) paste(x[1:2],collapse="-")))

fam <- metadata$family[idx]
pop <- metadata$PopulationName[idx]
sampname <- metadata$SampleName[idx]

match(as.character(pop),metadata2$pop.name)

table(table(as.character(fam)))

#pop <- sample(pop)
pop.mat <- outer(pop,pop, FUN = "==")
fam.mat <- outer(fam,fam, FUN = "==")

MRM(as.dist(as.matrix(dist01)[keep,keep]) ~ as.dist(pop.mat[keep,keep])
                       +as.dist(fam.mat[keep,keep]) 
          +   log1p(as.dist(as.matrix(geo.dist)[keep,keep])) )

# exclude pop == "LB-7"
keep <- pop != "LB-7"
fam.rep <- names(which (table(fam) > 1))
keep <- fam %in% fam.rep
MRM(as.dist(as.matrix(dist01)[keep,keep]) ~ 
#     as.dist(pop.mat[keep,keep])
#    as.dist(fam.mat[keep,keep]) 
    log1p(as.dist(as.matrix(geo.dist)[keep,keep])) )

plot( log1p(as.dist(as.matrix(geo.dist)[keep,keep])) , as.dist(as.matrix(dist01)[keep,keep]), main = "Genomic Isolation by Geographic Distance")

plot(density(dist01),xlim = c(0,1),main = "Distribution of pairwise genetic distance")
lines(density(dist01[as.dist(fam.mat)==1]),col="green")
lines(density(dist01[as.dist(pop.mat)==1]),col="blue")
legend(0,5,legend = c("all","family","population"),col = c("black","green","blue"), text.col = c("black","green","blue"),lty=1)

## load bioclim
library(raster)
library(rgdal)
xxx <- getData('worldclim',var='bio',res=2.5)

colors <- as.numeric(factor(as.character(pop)))

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
hc <- dendrapply(as.dendrogram(hc),function (n) attr(n, "nodePar") <- c(a$nodePar, lab.col=colors[n], pch=NA))
par(cex=0.3)
plot(hc)
