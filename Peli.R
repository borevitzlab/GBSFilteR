# R code for Population Genomics of Aus Pelargonia
# Justin Borevitz
# also used in Biol3157 practical 2014

# Genotyping By Sequencing

# HapMapFilter.R
# to run from your computer first
# download R from http://r-project.org

# use the Pelargonia example to start REQUIRED
snpFile <- "http://borevitzlab.anu.edu.au/~borevitz/gbs/HapMap.pel2bn3.hmn.txt"
## this file is a recall of only P. australae
snpFile <- "http://gduserv.anu.edu.au/~borevitz/tassel/pel1to3_paper/hapMap/HapMap.hmn.txt"

# if you are brave, try this Brachypodium dataset for ASSIGNMENT 2
#snpFile <- "http://borevitzlab.anu.edu.au/~borevitz/gbs/HapMap.bd2346.hmn.txt"

# look here for other pre-release data sets to play with OPTIONAL
# http://borevitzlab.anu.edu.au/~borevitz/gbs/

geno <- read.table(snpFile,header=T,na.strings=".")
# drop SNP info columns
g <- geno[,-(1:11)]

# this section matched plate and well with MetaData file
# renaming  for pel2bn3 data set only 
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
row.0 <- sub("0","",names.mat[5,])
tens <- grep("[[:upper:]]1",names.mat[5,])
row.0[tens] <- names.mat[5,tens]

#MetaData file
nameFile <- "http://borevitzlab.anu.edu.au/~borevitz/gbs/PeliSampleData-SpecimenWell.csv"
pel.names <-read.csv(nameFile)

index <- match(paste(names.mat[3,],row.0), paste(pel.names$Plate.Number,pel.names$Well.Number))
names.mat <- rbind(names.mat, 
	as.character(pel.names$Well.Number[index]),
	as.character(pel.names$Locality[index]))

# rename columns
names(g) <- paste(names.mat[1,],names.mat[7,])
#### ABOVE FOR PELI ONLY

# skip this step for now.
# drop Peli known contaminants identified previously
 contaminants <- names(g) %in%
    c("af.Bowk.3 na","af.Crisp.1 na","aus.D.2.34 Mt Chudalup","aus.AC.20.13 Rosedale","aus.I.18.06 Mt Stromlo") 
g <- g[,!contaminants]
dim(g)

## start saving output files
# pdf(file = "myNAMEassignment.pdf") # edit myNAME
# at the end of the script you need
# dev.off()

# Drop bad samples
samples<-apply(g,2,function(x) sum(!is.na(x)));
hist(samples,breaks = ncol(g),main = "Pel2bn3 samples with many SNPs",xlab = "SNPs")
# choose your own threshold here!
sample.snp.cutoff <- 1100
abline(v=sample.snp.cutoff,col = "red",lty=2)
table(samples >sample.snp.cutoff)
gg <- g[,samples >sample.snp.cutoff]
names.mat <- names.mat[,samples >sample.snp.cutoff]
table(is.na(gg))

## EXTRACT some info about the samples we will finally use from the MetaData
SampleName <- paste(pel.names$Complex,pel.names$Entity,pel.names$Plant.Number,sep=".")
use.index <- match(names(gg),paste(SampleName,pel.names$Locality))
table(pel.names$Complex[use.index])
table(pel.names$EntityCode[use.index])
table(pel.names$EntityCode[use.index],pel.names$Complex[use.index])

## trim low coverage SNPs
sites <- apply(gg,1, function(x) sum(!is.na(x)));
hist(sites,breaks=ncol(gg),main = "number of snps per sample", xlab = "samples")

# you pick another threshold
site.thresh.low <- 10 #Justin Borevitz
site.thresh.high <- 150 #Justin Borevitz

abline(v= c(site.thresh.low,site.thresh.high),lty=2)
site.index <- sites > site.thresh.low & sites < site.thresh.high
# list your output
table(site.index)
FALSE  TRUE 
23301  6230 

g.minor <- gg[site.index,]
dim(g.minor)

# what output did you get?
# [1] 6174  178
# rare variants
# Look for heterozygous rare variants
rare.var <- rowSums(g.minor,na.rm=T) == 1
table(rare.var)
# what output did you get?
#rare.var
# FALSE  TRUE 
# 6147    27 

# remove them if any
g.minor <- g.minor[!rare.var,]

# look for homozygous rare variants
rare.var <- rowSums(g.minor==2,na.rm=T) == 1
g.minor <- g.minor[!rare.var,]
# what output did you get?

rare.var <- rowSums(g.minor,na.rm=T) == 2
table(rare.var)
g.minor <- g.minor[!rare.var,]

# what output did you get?
#rare.var
#FALSE  TRUE 
# 5483    67 

image(as.matrix(g.minor),col = c("red","orange","blue"))

# save before removing paralogs with genotype ==1 color 
image(cbind(c(NA,0,1,2)),col = c("red","orange","blue"))
 # color coding, white is NA, red is 0, orange is 1, yellow is 2
paralog.g.minor <- g.minor

# paralogs/duplicate markers 
het.count <- rowSums(g.minor==1,na.rm=T)
hist(het.count,breaks=ncol(g.minor))
paraThresh <- 90
abline(v = paraThresh, lty = 2)
table(het.count < paraThresh)

# what output did you get?
FALSE  TRUE 
   24  5526 

# remove paralogs
g.minor <- g.minor[het.count < paraThresh,]

# notice the paralogs are removed.
het.count <- rowSums(g.minor==1,na.rm=T)
hist(het.count,breaks=ncol(g.minor))
# close raw plot file
#dev.off()

table(is.na(g.minor))

# how many total calls?
# FALSE   TRUE 
#238660 744968 


# same images as .jpg to make them smaller
png(width = 1024, height = 768, file = "PelSNP.png")
image(as.matrix(paralog.g.minor))
dev.off()

#don't need to save out raw snp table
write.csv(g.minor, file = "Pel2bn3.5527.snp.181samp.csv")

# dendrogram
corDist <- as.dist(1-cor(g.minor,use="pairwise.complete.obs")^2)
hc <- hclust(corDist)

d.dist <- dist(t(g.minor))
hc <- hclust(d.dist)

# Loop to find and drop pairs of samples that are returning missing values
# Probably not necessary if “af.Bowk.3 na” was dropped as a contaminant in Pel data
# Use with caution.
# Author: Ben Bai
# g.minor.filtered <- g.minor
#
# repeat {
# 	samples.to.drop <- unique(rownames(which(is.na(as.matrix(d.dist)), arr.ind=T)))
# 	print("Dropping samples causing missing distance matrix values:")
# 	print(samples.to.drop)
# 	if (length(samples.to.drop) == 0) {
#     		break
# 	}
# 	g.minor.filtered <- g.minor.filtered[,!(names(g.minor) %in% samples.to.drop)]
# 	d.dist <- dist(t(g.minor.filtered))
# }
#
# d.dist.filtered = dist(t(g.minor.filtered))
# hc <- hclust(d.dist.filtered)

# rename with your SNPnumber and student name
#pdf(file = "pel2b3n.mySNPs.mySamples.pdf",height = 8, width = 11)
par(cex=0.6)# adjust text size
plot(hc)
#dev.off()

# describe a finding in the .pdf you make below, for example which sites names cluster and which don’t

## look far below for the last PCoA plots.
# examples 
#We saw that lit, rod, and aus were the major deep splits in the tree, less defined clustering was seen at the lower levels. Example, ausAC and ausAI are more mixed. The other points from the African samples were seen to cluster around the aus samples. #150 Kira Anderson

#Clusters: lit.Ha, lit. L, rod.S/rod.R/aus.D/aus.AC, af.Gross/af.lon/af.ren, af.Killing, aus.H/aus.I, aus.IA/aus.AC, aus.AI/aus.AC/aus.I, aus.D, af.Crisp/aus.AC/aus.AI/rod.AI, these clustering shows that the splits are grouped accordingly to the family, which one or two species crossing over to another cluster. #90 Nay Chi Khin


# Principal Coordinates
test <- t(as.matrix(g.minor))
test[test==0] <- "c"
test[test==1] <- NA
test[test==2] <- "t"

#install.packages("ape") # if package was not installed
library(ape)
hap.genotypes.bin <- as.DNAbin(test)
hap.gene.dist <- dist.dna(hap.genotypes.bin, model = "raw" , pairwise.deletion=T)
#image(as.matrix(hap.gene.dist)) #not that important to include
par(cex=0.6)# adjust text size
plot(hclust(hap.gene.dist))

# multidimensional scaling.  Principal Coordinate analysis
pcoa <- cmdscale(hap.gene.dist, k = 8)
pcoa_df <- data.frame(pcoa)
tot.var <- sum(sapply(pcoa_df, var))
per.var <- round((sapply(pcoa_df, var)/tot.var)*100,2)

# try to color names by locality but not
pcoa.color <- rep(1,ncol(g.minor)) #
pcoa.color[grep("lit",names(g.minor))] <- 3
pcoa.color[grep("aus",names(g.minor))] <- 4
pcoa.color[grep("rod",names(g.minor))] <- 2

# for Brachy bd2346
#pcoa.color[grep("Bdis",names(g.minor))] <- 2
#pcoa.color[grep("BdTR",names(g.minor))] <- 3
# PC1 vs PC2

plot(pcoa[,2]~pcoa[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Pel2bn3 PCoA" , col = pcoa.color )
text(pcoa[,2] ~ pcoa[,1], 
                            labels =colnames(g.minor) ,
                            cex=0.5,col = pcoa.color )
legend(min(pcoa[,1]), max(pcoa[,2]),c("lit","aus","rod","other"),col= c(3,4,2,1),lty=1)

dev.off() # close off the .pdf to finish assignment

# PC2 vs PC3
#pdf(file = "pel2bn3myNamePCoA.pdf",paper = "a4r",width = 11, height = 15)

plot(pcoa[,3]~pcoa[,2], xlab = paste("PC2 var", per.var[2]),
                            ylab = paste("PC3 var", per.var[3]),
                            main = "Pel2bn3 PCoA" , col = pcoa.color )
text(pcoa[,3] ~ pcoa[,2], 
                            labels =colnames(g.minor) ,
                            cex=0.5,col = pcoa.color )
legend(min(pcoa[,2]), max(pcoa[,3]),c("lit","aus","rod","other"),col= c(3,4,2,1),lty=1)
#legend(min(pcoa[,1]), max(pcoa[,2]),c("Bdis","BdTR","other"),col= c(2,3,1),lty=1)
dev.off() 


