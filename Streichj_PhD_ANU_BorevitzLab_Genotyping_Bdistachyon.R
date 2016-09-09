################################################################################################
## <---------------------------- Window Width for Script -----------------------------------> ##
########################################## Citation ############################################
# Genotyping Brachypodium distachyon: PhD project of Jared Streich at the Borevitz lab at ANU
#  Version D.2.0
#  Authors: Jared Streich, Justin Borevitz, Kevin Murray, Steve Eichten, 
#           Aaron Chuah, Jason Bragg

########################################### CONTENTS ###########################################
# Packages to Install ----------------------------------------------- 24
# Load Data - Read in Data Frames ----------------------------------- 113
# GENETIC - hapmap anaylysis ---------------------------------------- 235

################################## Set Directory ###############################################
setwd("~/Desktop/R_files")

################################################################################################
################################# Call Packages for Script #####################################
################################################################################################

library(phangorn)
require(ggplot2)
library(ape)
library(ggdendro)
library(graphics)
library(cluster) 
library(dynamicTreeCut)
library(reshape2)
library(stats)
library(plyr)
library(RFLPtools)
library(adephylo)
library(dendextend)
library(mclust)

################################################################################################
############################## Install packages if needed ######################################
################################################################################################

install.packages("mclust")
install.packages("car")
install.packages("adephylo")
install.packages("dynamicTreeCut")
install.packages('ggplot2')
install.packages('ggdendro')
install.packages("phangorn")
install.packages('RFLPtools')
install.packages("e1071")
install.packages("apcluster")


########### Clear all remaining variables in R, 
########### CAREFULL... you don't ruin another script running!!! 
rm(list = ls())

________________________________________________________________________________________________
________________________________________________________________________________________________

# LLL      OOOOOO  AAAAAA    DDDDDD     DDDDDDD  AAAAAA    TTTTTTTTTTT AAAAAA    
# LLL     OOO  OOO AAA AAA   DDDDDDDD   DDDDDDDD AAA AAA   TTTTTTTTTTT AAA AAA   
# LLL     OOO  OOO AAA  AAA  DDD  DDD   DDD  DDD AAA  AAA      TTT     AAA  AAA  
# LLL     OOO  OOO AAAAAAAAA DDD  DDD   DDD  DDD AAAAAAAAA     TTT     AAAAAAAAA 
# LLL     OOO  OOO AAAAAAAAA DDD  DDD   DDD  DDD AAAAAAAAA     TTT     AAAAAAAAA 
# LLLLLLL OOO  OOO AAA   AAA DDDDDDDD   DDDDDDDD AAA   AAA     TTT     AAA   AAA 
# LLLLLLL  OOOOOO  AAA   AAA DDDDDDD    DDDDDDD  AAA   AAA     TTT     AAA   AAA 
________________________________________________________________________________________________
________________________________________________________________________________________________

#######################################################################################
################################## Read Hapmap Files ##################################
#######################################################################################

# All Genomic Chromosomes
geno <- read.table("kmeansDistachyon.txt",header=T,na.strings=".")
dim(geno)
head(geno)

################# Strip marker info from Header
# Remove Variant Data Genomic Data
g <- geno[,-(1:11)] 
# Remove Duplicate data from memory
rm(geno)

# All Genomic Data Format Names For Data Sheet Matching
names(g)[1] <- "ref-Dgenome_12_12_12_12"
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
colnames(g) <- names.mat[1,]


#######################################################################################
############################### Read Meta Data Files ##################################
#######################################################################################

################### Read in Meta Data 
Meta.Data <- read.table("Brachy_META_NO_NAs.txt", header=T, sep="	")
rownames(Meta.Data) <- Meta.Data[,1]
head(Meta.Data)
dim(Meta.Data)

#######################################################################################
########################### Remove Suspect and bad sampels ############################
#######################################################################################
# Many of these are known polyploids or samples that had
# intermediate proportions of B. stacei markers and were 
# flagged and removed, or artifacts from a known detectable 
# glitch in AXE demultiplexing output, A few are also from
# a plate that had errors and barcode problems and should
# be resequenced.

colnames(g)
dim(g)
strt.g <- ncol(g)

g <- g[,-522] # Luc9
g <- g[,-520] # Luc3
g <- g[,-518] # Luc21
g <- g[,-499] # KJP.7
g <- g[,-498] # KJP.3
g <- g[,-497] # KJP.2
g <- g[,-496] # KEI.5
g <- g[,-495] # KEI.1
g <- g[,-494] # KAH.3
g <- g[,-491] # Isk.p8
g <- g[,-490] # Isk.p7
g <- g[,-489] # Isk.p6
g <- g[,-488] # Isk.p4
g <- g[,-487] # Isk.p2
g <- g[,-479] # Gaz6
g <- g[,-469] # Gal6.2
g <- g[,-439] # Bel15.2
g <- g[,-437] # Bei3.2
g <- g[,-436] # Bdis403.1
g <- g[,-435] # Bdis402.1
g <- g[,-364] # Bdis28.19.1
g <- g[,-311] # BdTR11e
g <- g[,-216] # BdTR9f
g <- g[,-153] # BdTR2a
g <- g[,-146] # BdTR1g
g <- g[,-141] # BdTR1c
g <- g[,-137] # BDTR.13o
g <- g[,-134] # BdTR13k
g <- g[,-133] # BdTR13d
g <- g[,-115] # BdTR11e
g <- g[,-112] # BdTR11b
g <- g[,-109] # BdTR10n
g <- g[,-96]  # BdTR.8c
g <- g[,-93]  # Bd21
g <- g[,-92]  # Bd21.1
g <- g[,-91]  # Bd21.1.2 
g <- g[,-90]  # Bd21.3 
g <- g[,-64]  # BDTR.13o
g <- g[,-62]  # BD3.1
g <- g[,-59]  # BD2.3
g <- g[,-57]  # BD19.6
g <- g[,-56]  # BD13.8
g <- g[,-14]  # ABR7.R3
g <- g[,-13]  # ABR6.R1
g <- g[,-11]  # ABR5.R4
g <- g[,-5]   # X51SPA.Nor.s6E.c3.R4
g <- g[,-4]   # X18Geo32i2.C4.R4

rmvd.g <- ncol(g)
fnl <- strt.g-rmvd.g

# Total samples removed so far
fnl

________________________________________________________________________________________________
________________________________________________________________________________________________

#   GGGGGGG  EEEEEEEE NNNN   NNN EEEEEEEEE TTTTTTTTTTTT IIIIIIIIII   CCCCCCCC
#  GGG   GGG EEE      NNNNN  NNN EEE       TT  TTTT  TT II IIII II  CCC   CCCC
#  GGG       EEEEEEEE NNNNNN NNN EEEEEEEEE     TTTT        IIII     CCC   
#  GGG  GGGG EEEEEEEE NNN NNNNNN EEEEEEEEE     TTTT        IIII     CCC
#  GGG  G GG EEE      NNN  NNNNN EEE           TTTT        IIII     CCC
#  GGG    GG EEEEEEEE NNN   NNNN EEEEEEEEE     TTTT     II IIII II  CCC   CCCC
#   GGGGGGG  EEEEEEEE NNN    NNN EEEEEEEEE     TTTT     IIIIIIIIII   CCCCCCCC

________________________________________________________________________________________________
________________________________________________________________________________________________

################################################################################################
############### This section is for analysing hapmap data and creating figures #################
################################################################################################

################################################################################################
######################### First Glance at Raw Hapmap Data ######################################
################################################################################################

############ Image raw data
png(file="Kmeans_Distachyon_Raw_SNP_DataImage.png",width = 3072, height = 3*768)
image(1:nrow(g),1:ncol(g),as.matrix(g),xlab = "SNPs", ylab = "samples", 
       main = paste(c("SNPs","samples"),dim(g)) )
dev.off()

################################################################################################
############################### Start Filtering ################################################
################################################################################################

# minimum markers for sample use
mm <- 10000


pdf(file = "Distachyon_Histogram_SampleCalls.pdf",width = 30, height = 30)
######## Identify the number of genotype calls by marker
samplecalls<-apply(g,2, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls,breaks=ncol(g))
samp.thres <- mm
abline(v=samp.thres,lty=2)
dev.off()
mean(samplecalls)



# minimum shared alleles
ma <- max(snpsPerSample)*0.55
ma


pdf(file = "Distachyon_genome_Histogram_SNPsPerSample.pdf",width = 30, height = 30)
################# Test a couple thresholds
table(samplecalls > samp.thres)
g.snp <- g[,samplecalls > samp.thres]
############## Look at counts across samples
snpsPerSample<-apply(g.snp,1,function(x) sum(!is.na(x)));
hist(snpsPerSample,breaks = 100)
############## Test some sample count thresholds
snpCallThresh <- ma
abline(v = snpCallThresh,lty= 2)
table(snpsPerSample > snpCallThresh)
g.minor <- g.snp[snpsPerSample > snpCallThresh,]
dev.off()


######################################### Rare varients ########################################
#### 1 het
rare.var <- rowSums(g.minor,na.rm=T) == 1
table(rare.var)
g.minor <- g.minor[!rare.var,]

##### 1 hom
rare.var <- (rowSums(g.minor==2,na.rm=T) == 1) & (rowSums(g.minor,na.rm=T) == 2)
table(rare.var)
g.minor <- g.minor[!rare.var,]

############### Paralogs?
#image(g.minor==1)

pdf(file = "Distachyon_Histogram_HetCounts.pdf",width = 30, height = 30)
######## paralogs/duplicate markers 
het.count <- rowSums(g.minor==1,na.rm=T)
hist(het.count,breaks=ncol(g.minor))
paraThresh <- 5
abline(v = paraThresh, lty = 2)
table(het.count < paraThresh)
dev.off()


g.minor <- g.minor[het.count < paraThresh,]
########## How many total calls? ################
table(is.na(g.minor))
dim(g.minor)



########### minor allele frequency distribution
hist(rowSums(g.minor,na.rm=T)/2, 100)

########### Transform Matrix
test <- t(as.matrix(g.minor))
test[test==0] <- "c"
test[test==1] <- NA
test[test==2] <- "t"

# Create Priciple Coordinate Vectors
hap.genotypes.bin <- as.DNAbin(test)
hap.gene.dist <- dist.dna(hap.genotypes.bin, model = "N" , pairwise.deletion=T)
swil.1 <- cmdscale(hap.gene.dist, k = 19)
swil.1_df <- data.frame(swil.1)
tot.var <- sum(sapply(swil.1_df, var))
per.var <- round((sapply(swil.1_df, var)/tot.var)*100,2)

# Create Pairwise Distance for Dendrograms
snp.1cor2.dist <- snp.1cor2.dist <- as.dist(1-cor(g.minor,use = "pairwise.complete.obs")^2)
dim(swil.1)


###############################################################################################
###################### Format Meta Data to Genetic info for plotting ##########################
###############################################################################################

# Match Names
MD.g <- match(rownames(Meta.Data), rownames(swil.1))
MD.Dist <- cbind(Meta.Data, MD.g)

# Check Dimensions and Names
head(MD.Dist)
dim(MD.Dist)

# Filter Meta Data by Samples that didn't match
MD.Dist <- MD.Dist[complete.cases(MD.Dist$MD.g),]

# Remove duplicate information form memory
rm(MD.g)
rm(Meta.Data)

# Check the dimensions of new Meta Data Matrix
head(MD.Dist)
dim(MD.Dist)
rownames(MD.Dist)
match(rownames(swil.1), rownames(MD.Dist))
dim(swil.1)
MD.Dist <- MD.Dist[order(rownames(MD.Dist)),]
swil.1 <- swil.1[order(rownames(swil.1)),]

################################################################################################
############################## Create Visuals Post Filtering ###################################
################################################################################################




###################### Principle Coordinate by Plate ##################
pdf(file = "Distachyon_Plate_PCoA.pdf",width = 30, height = 30)
palette(rainbow(max(MD.Dist$Plate)))
plot(swil.1[,2], swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachypodium distachyon PCoA by Plate", 
                            col=MD.Dist$Plate, cex = 2, pch = 16)
text(swil.1[,2], swil.1[,1], labels = names(g), cex=0.5,col="Black")
legend(-1200,400,paste("Plate",c(1:32)),
							col = c(1:32),lty=1, cex = 0.7)
dev.off()



###################### Principle Coordinate by Plate ##################
pdf(file = "Distachyon_Continent_PCoA.pdf",width = 30, height = 30)
palette(rainbow(4))
plot(swil.1[,2],swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachypodium distachyon PCoA by Plate", 
                            col=MD.Dist[,7], cex = 5, pch = 16)
text(swil.1[,2], (swil.1[,1]), labels = colnames(g), cex=0.5,col="Black")
legend(-1200,300,paste("Continent",c(1:4)),
							col = c(1:4),lty=1, cex = 3)
dev.off()


######################### PDF of Kinship, Heatmap ########################
pdf(file = "Distachyon_Genetic_Kinship.pdf", width = 30, height = 30)
image(1:nrow(test),1:nrow(test),as.matrix(dist.dist))
dev.off()


#################### Count Groups by Unrooted Tree ############################
CG <- hclust(dist(snp.1cor2.dist))
plot(as.phylo(CG), type = "unrooted")

#  Define Number of Major Groups
G.groups = 101
G.height = 1
G.mclust = 16
palette(rainbow(G.groups))


############ MClust for Genotping, or at least major families #################
Gene.mclust <- Mclust(as.matrix(swil.1))
G.best <- dim(Gene.mclust$z)[2]
cat("model-based optimal number of clusters:", G.best, "\n")
plot(Gene.mclust)

dim(g.minor)
################################ Calinksi Best K ##################################

half.samps <- (dim(g.minor)/2)
Gene.fit <- cascadeKM(scale(g.minor, center = TRUE,  scale = TRUE), 2, half.samps, iter = 1000)
plot(Gene.fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(Gene.fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")




########################### Test Dendrogram no hang ###########################
pdf(file = "Distachyon_101UniqueGroups_256samps_min33=8kSNPs_DendroNoHang.pdf",width = 50, height = 20)
plot((CG), cex=0.5, hang = -1)
rect.hclust(CG, k = G.groups, border = "Red")
dev.off()

######################### PDF of Dendrogram no hang By Groups #######################
pdf(file = "Distachyon_height-1_UniqueGroups_DendroNoHang.pdf",width = 50, height = 20)
plot((CG), cex=0.3)
rect.hclust(CG, h=G.height, border = "Red")
dev.off()

###################### Output names of Dendrogram by Groups ######################
write.hclust(CG, file="Distachyon_101Unique_Groups.txt", prefix = "BdDendro", 
								h = NULL, k = G.groups, append = FALSE, dec = ",")



# Read and Match Names by Group Cluster
Clust.101.Groups <- read.delim(file="Distachyon_101Unique_Groups.txt", header = T)
head(Clust.101.Groups)
rownames(Clust.101.Groups) <- Clust.101.Groups[,1]
match(rownames(MD.Dist), rownames(Clust.101.Groups))
head(Clust.101.Groups$Cluster)



# Output at four Groups to test hypothetical K value for pie charts on maps, at k = 4
write.hclust(CG, file="Distachyon_14Unique_Groups.txt", prefix = "BdDendro", 
								h = NULL, k = 14, append = FALSE, dec = ",")



# Read and Match Names by Group Cluster
Clust.14.Groups <- read.delim(file="Distachyon_14Unique_Groups.txt", header = T)
head(Clust.14.Groups)
rownames(Clust.14.Groups) <- Clust.Groups[,1]
match(rownames(MD.Dist), rownames(Clust.14.Groups))
match(rownames(swil.1), rownames(Clust.14.Groups))
head(Clust.14.Groups$Cluster)

######################## TEst Plot for Pie groups at 14 colors
palette(rainbow(14))
plot(swil.1[,2], swil.1[,1], col = Clust.14.Groups[,2], cex = log(swil.1[,3])+1, lwd = 3)


################################################################################################
############################### For height specific genotype calling ###########################
################################################################################################

######################### PDF of Dendrogram no hang By Height #######################
pdf(file = "Distachyon_479Samps_131Groups_DendroNoHang.pdf",width = 50, height = 20)
plot((CG), cex=0.3)
rect.hclust(CG, h = 1, border = "Red")
dev.off()

# Read and Match Names of Groups by Height Cluster
write.hclust(CG, file="Distachyon_Groups_h=1.txt", prefix = "BdDendro", 
		h = 1, append = FALSE, dec = ",")

# Read and Match Names by Group Cluster by Height
Clust.Height.1 <- read.delim(file="Distachyon_Groups_h=1.txt", header = T)
head(Clust.Height.1)
rownames(Clust.Height.1) <- Clust.Height.1[,1]
match(rownames(MD.Dist), row.names(Clust.Height.1))

# Number of Unique Groups at height 1
max(Clust.Height.1$Cluster)


# Read and Match Names of Groups by Height Cluster
write.hclust(CG, file="Distachyon_Groups_h=2.5.txt", prefix = "BdDendro", 
		h = 2.5, append = FALSE, dec = ",")

# Read and Match Names by Group Cluster by Height
Clust.Height.2.5 <- read.delim(file="Distachyon_Groups_h=2.5.txt", header = T)
head(Clust.Height.2.5)
row.names(Clust.Height.2.5) <- Clust.Height.2.5[,1]
match(rownames(MD.Dist), row.names(Clust.Height.2.5))

# Number of Unique Groups at height 2.5
max(Clust.Height.2.5$Cluster)




# Read and Match Names of Groups by Height Cluster
write.hclust(CG, file="Distachyon_3Groups.txt", prefix = "BdDendro", 
		k = 3, append = FALSE, dec = ",")

# Read and Match Names by Group Cluster by Height
Clust.3.groups <- read.delim(file="Distachyon_3Groups.txt", header = T)
head(Clust.3.groups)
row.names(Clust.3.groups) <- Clust.3.groups[,1]
match(rownames(MD.Dist), row.names(Clust.3.groups))

# Number of Unique Groups at height 2.5
max(Clust.3.groups$Cluster)




pdf(file = "Distachyon_G.groups_PCoA.pdf",width = 30, height = 30)
plot(swil.1[,2]~swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachy PCoA", cex = 1, pch = Clust.Height.1$Cluster)
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels = names.mat ,
                            cex=0.5,col="Black")
dev.off()




################################################################################################
############################## Create Structure Input File #####################################
################################################################################################

# Transform SNP Table
Struct.Dist <- t(as.matrix(g.minor))
Struct.Dist.Names <- rownames(Struct.Dist)
head(Struct.Dist)
dim(Struct.Dist)

Struct.Dist[is.na(Struct.Dist)] <- -9

# Re-Attach Row Names
row.names(Struct.Dist) <- Struct.Dist.Names

# Output Structure File
write.table(Struct.Dist, "Struct.Dist.txt", sep="	", row.names = T, col.names = T)
dim(Struct.Dist)
head(Struct.Dist)
Dist.sum <- Struct.Dist

Dist.sum[Dist.sum == -9] <- 0
head(Dist.sum)
Dist.col.sum <- colSums(Dist.sum)


pdf(file = "Distachyon_Marker_Sums.pdf",width = 100, height = 15)
plot(Dist.col.sum, cex = 0.2, type = "b", pch = 16)
dev.off()

