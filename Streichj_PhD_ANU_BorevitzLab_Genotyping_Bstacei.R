################################################################################################
## <---------------------------- Window Width for Script -----------------------------------> ##
########################################## Citation ############################################
# Genotyping Brachypodium stacei: the PhD project of Jared Streich at the Borevitz lab at ANU
#  Version S.2.0
#  Authors: Jared Streich, Justin Borevitz, Kevin Murray, Steve Eichten, 
#           Aaron Chuah, Jason Bragg

########################################### CONTENTS ###########################################
# Packages to Install ----------------------------------------------- 
# Load Data - Read in Data Frames ----------------------------------- 
# GENETIC - hapmap anaylysis ----------------------------------------

################################## Set Directory ###############################################
setwd("~/Desktop/")

################################################################################################
################################# Call Packages for Script #####################################
################################################################################################

library(phangorn)
require(ggplot2)
library(ape)
library(ggdendro)
library(graphics)
library(cluster) 
require(grDevices)
library(dynamicTreeCut)
library(reshape2)
library(stats)
library(plyr)
library(RFLPtools)
library(adephylo)
library(mclust)


################################################################################################
############################## Install packages if needed ######################################
################################################################################################

install.packages("mclust")
install.packages("pegas")
install.packages("car")
install.packages("maps")
install.packages("adephylo")
install.packages("synbreed")
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
geno <- read.table("kmeansStacII.txt",header=T,na.strings=".")
dim(geno)
head(geno)

################# Strip marker info from Header
# Remove Variant Data Genomic Data
g <- geno[,-(1:11)] 
# Remove Duplicate data from memory
rm(geno)

# All Genomic Data Format Names For Data Sheet Matching
names(g)[1] <- "ref-Sgenome_12_12_12_12"
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
png(file="Stacei_Raw_SNP_DataImage.png",width = 3072, height = 3*768)
image(1:nrow(g),1:ncol(g),as.matrix(g),xlab = "SNPs", ylab = "samples", 
       main = paste(c("SNPs","samples"),dim(g)) )
dev.off()

################################################################################################
############################### Start Filtering ################################################
################################################################################################


pdf(file = "stacei_Histogram_SampleCalls.pdf",width = 30, height = 30)
######## Identify the number of genotype calls by marker
samplecalls<-apply(g,2, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls,breaks=ncol(g))
samp.thres <- 6000
abline(v=samp.thres,lty=2)
dev.off()

mean(samplecalls)

pdf(file = "Stacei_genome_Histogram_SNPsPerSample.pdf",width = 30, height = 30)
################# Test a couple thresholds
table(samplecalls > samp.thres)
g.snp <- g[,samplecalls > samp.thres]
############## Look at counts across samples
snpsPerSample<-apply(g.snp,1,function(x) sum(!is.na(x)));
hist(snpsPerSample,breaks = 100)
############## Test some sample count thresholds
snpCallThresh <- 17
abline(v = snpCallThresh,lty= 2)
table(snpsPerSample > snpCallThresh)
g.minor <- g.snp[snpsPerSample > snpCallThresh,]
dev.off()

max(snpsPerSample)*0.33

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

pdf(file = "Stacei_Histogram_HetCounts.pdf",width = 30, height = 30)
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
g2.minor <- g.minor[23185:24185,]


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


################################################################################################
############################## Create Visuals Post Filtering ###################################
################################################################################################


###################### Principle Coordinate by Plate ##################
pdf(file = "Stacei_Plate_PCoA.pdf",width = 10, height = 10)
palette(rainbow(max(MD.Dist$Plate)))
plot(swil.1[,2], swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachypodium stacei PCoA by Plate", 
                            col=MD.Dist$Plate, cex = 2, pch = 16)
text(swil.1[,2], swil.1[,1], labels = names(g), cex=0.5,col="Black")
legend(-1200,400,paste("Plate",c(1:32)),
							col = c(1:32),lty=1, cex = 0.7)
dev.off()





###################### Principle Coordinate by Plate ##################
pdf(file = "stacei_Continent_PCoA.pdf",width = 10, height = 10)
palette(rainbow(4))
plot(swil.1[,2],swil.1[,1], xlab = paste("PC1 var", per.var[1]),
                            ylab = paste("PC2 var", per.var[2]),
                            main = "Brachypodium stacei PCoA by Plate", 
                            col=MD.Dist$Continent, cex = 5, pch = 16)
text(swil.1[,2], (swil.1[,1]), labels = colnames(g), cex=0.5,col="Black")
legend(-1200,300,paste("Continent",c(1:4)),
							col = c(1:4),lty=1, cex = 3)
dev.off()


######################### PDF of Kinship, Heatmap ########################
pdf(file = "stacei_Genetic_Kinship.pdf", width = 30, height = 30)
image(1:nrow(test),1:nrow(test),as.matrix(dist.dist))
dev.off()


#################### Count Groups by Unrooted Tree ############################
CG <- hclust(dist(snp.1cor2.dist))
plot(as.phylo(CG), type = "unrooted")


#  Define Number of Major Groups
G.groups = 8
G.height = 1
G.mclust = 16
palette(rainbow(G.groups))


############ MClust for Genotping, or at least major families #################
Gene.mclust <- Mclust(as.matrix(swil.1))
G.best <- dim(Gene.mclust$z)[2]
cat("model-based optimal number of clusters:", G.best, "\n")
plot(Gene.mclust)


################################ Calinksi Best K ##################################
Gene.fit <- cascadeKM(scale(swil.1, center = TRUE,  scale = TRUE), 2, 51, iter = 1000)
plot(Gene.fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(Gene.fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")




########################### Test Dendrogram no hang ###########################
pdf(file = "stacei_13UniqueGroups_no-hang_DendroNoHang.pdf",width = 50, height = 20)
plot((CG), cex=0.5, hang = -1)
rect.hclust(CG, k = G.groups, border = "Red")
dev.off()

######################### PDF of Dendrogram no hang By Groups #######################
pdf(file = "stacei_8UniqueGroups_DendroNoHang.pdf",width = 50, height = 20)
plot((CG), cex=0.3)
rect.hclust(CG, G.groups, border = "Red")
dev.off()

###################### Output names of Dendrogram by Groups ######################
write.hclust(CG, file="Stacei_13Unique_Groups.txt", prefix = "BdDendro", 
								h = NULL, k = G.groups, append = FALSE, dec = ",")



# Read and Match Names by Group Cluster
Clust.Groups <- read.delim(file="Stacei_13Unique_Groups.txt", header = T)
head(Clust.Groups)
rownames(Clust.Groups) <- Clust.Groups[,1]
match(rownames(MD.Dist), rownames(Clust.Groups))
head(Clust.Groups$Cluster)



# Output at four Groups to test hypothetical K value for pie charts on maps, at k = 4
write.hclust(CG, file="stacei_4Unique_Groups.txt", prefix = "BdDendro", 
								h = NULL, k = 16, append = FALSE, dec = ",")



# Read and Match Names by Group Cluster
Clust.4.Groups <- read.delim(file="stacei_4Unique_Groups.txt", header = T)
head(Clust.Groups)
rownames(Clust.Groups) <- Clust.Groups[,1]
match(rownames(MD.Dist), rownames(Clust.Groups))
head(Clust.Groups$Cluster)



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
write.table(Struct.Dist, "Struct.Stacei.txt", sep="	", row.names = T, col.names = T)
dim(Struct.Dist)
head(Struct.Dist)
Dist.sum <- Struct.Dist

Dist.sum[Dist.sum == -9] <- 0
head(Dist.sum)
Dist.col.sum <- colSums(Dist.sum)


pdf(file = "Stacei_Marker_Sums.pdf",width = 100, height = 15)
plot(Dist.col.sum, cex = 0.2, type = "b", pch = 16)
dev.off()



