################################################################################################
## <---------------------------- Window Width for Script -----------------------------------> ##

## Script Starts at line 80

########################################## Citation ############################################
#  Genotyping Brachypodium hybridum: the PhD project of Jared Streich at the Borevitz lab at ANU
#  Version H.2.0
#  Authors: Jared Streich, Justin Borevitz, Kevin Murray, Steve Eichten, 
#           Aaron Chuah, Jason Bragg
#
########################################### CONTENTS ###########################################
# Packages to Install
# Load Data - Assign text files to variables
# GENETIC - hapmap anaylysis

################################################################################################
################################# Call Packages as Needed ######################################
################################################################################################
library(phangorn)
require(ggplot2)
require(ape)
library(ggdendro)
require(graphics)
library(cluster) 
require(grDevices)
library(dynamicTreeCut)
library(reshape2)
require(stats)
library(plyr)
library(RFLPtools)
library(adephylo)
library(e1071)
library(mclust)

################################################################################################
############################## Install packages as needed ######################################
################################################################################################

install.packages("mclust")
install.packages("adephylo")
install.packages("synbreed", dep = TRUE)
install.packages("dynamicTreeCut")
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
install.packages('ggplot2')
install.packages('ggdendro')
install.packages("phangorn")
install.packages('RFLPtools')
install.packages('plotrix')
install.packages("e1071")


########### Clear all remaining variables in R, 
########### CAREFULL... you don't ruin another script running!!! 
rm(list = ls())


################################################################################################
############################## Change Working Directory ########################################
################################################################################################

# You will need to edit this to your directory where you have the text files.
setwd("/Users/u5212257/Desktop/R_files")


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

############# Read in hapmap files

# All Genomic Chromosomes
H.geno <- read.table("kmeansHybridumH.txt",header=T,na.strings=".")
dim(H.geno)

#######################################################################################
############################# Per Species/SubGenome ###################################
#######################################################################################


##########Load Sub-Genome Calls
S.geno <- read.table("kmeansHybridumS.txt",header=T,na.strings=".")
D.geno <- read.table("kmeansHybridumD.txt",header=T,na.strings=".")

################ Check Dimensions and Headers
dim(H.geno)
dim(S.geno)
dim(D.geno)
head(H.geno)
head(S.geno)
head(D.geno)

# Set SNP Names / RowNames
rownames(H.geno) <- H.geno[,1]
rownames(S.geno) <- S.geno[,1]
rownames(D.geno) <- D.geno[,1]

# Species/Sub-genomic Data
H <- H.geno[,-(1:11)] 
S <- S.geno[,-(1:11)] 
D <- D.geno[,-(1:11)] 
rownames(S)

dim(H)
dim(S)
dim(D)
rm(H.geno)
rm(S.geno)
rm(D.geno)


# Number of Variants in B. stacei genome
Svar <- nrow(S)
Dvar <- nrow(D)
Hvar <- nrow(H)


# Format Names For Data Sheet Matching
names(H)[1] <- "ref-Hgenome_12_12_12_12"
names.mat.H <- matrix(unlist(strsplit(names(H),split="_")),nrow=5)
names(H) <- names.mat.H[1,]

names(S)[1] <- "ref-genome_12_12_12_12"
names.mat.S <- matrix(unlist(strsplit(names(S),split="_")),nrow=5)
names(S) <- names.mat.S[1,]

# D Genome Names Format
names(D)[1] <- "ref-genome_12_12_12_12"
names.mat.D <- matrix(unlist(strsplit(names(D),split="_")),nrow=5)
names(D) <- names.mat.D[1,]


#######################################################################################
############################### Read Meta Data Files ##################################
#######################################################################################

################### Read in Meta Data 
Meta.Data <- read.table("Brachy_META_NO_NAs.txt", header=T, sep="	")
rownames(Meta.Data) <- Meta.Data[,1]
dim(Meta.Data)
head(Meta.Data)

# Match Names
MD.H <- match(rownames(Meta.Data), colnames(H), incomparables=NULL)
MD.S <- match(rownames(Meta.Data), colnames(S), incomparables=NULL)
MD.D <- match(rownames(Meta.Data), colnames(D), incomparables=NULL)

# Combine H, S, and D Meta Data
MD.HSD <- cbind(Meta.Data, MD.H, MD.S, MD.D)

# Combine S and D Meta Data
MD.SD <- cbind(Meta.Data, MD.S, MD.D)

# Check for names that don't match and edit if necessary
MD.namecheck <- match(names(S), rownames(Meta.Data))
MD.namecheck <- cbind(names(S), MD.namecheck)

# Check Dimensions and Names
head(MD.SD)
dim(MD.SD)

# Filter Meta Data by Samples that didn't match
MD.SD <- MD.SD[complete.cases(MD.SD$MD.S),]

# Remove duplicate information form memory
rm(MD.H)
rm(MD.S)
rm(MD.D)
rm(Meta.Data)

# Check the dimensions of new Meta Data Matrix
head(MD.HSD)
dim(MD.SD)





#######################################################################################
##################################### Species Check ###################################
#######################################################################################


# Match Names of S and D Genomic Data
match(names(S), names(D), nomatch = NA, incomparables = NULL)


######################### Convert calls to Presence-Absence of Markers
# Convert S genome markers to Presence-Absence
S.binary.step1 <- apply(S, 2, function(y) as.numeric(gsub(0, 1, y)))
S.binary.step2 <- apply(S.binary.step1, 2, function(w) as.numeric(gsub(2, 1, w)))

# Convert D genome markers to Presence-Absence
D.binary.step1 <- apply(D, 2, function(z) as.numeric(gsub(0, 1, z)))
D.binary.step2 <- apply(D.binary.step1, 2, function(v) as.numeric(gsub(2, 1, v)))

# Sum Columns
S.sum <- colSums(S.binary.step2, na.rm = T)
D.sum <- colSums(D.binary.step2, na.rm = T)

######################## Confirm names are attached to columns
names(S.sum)=names(S)
names(D.sum)=names(D)

# Plot B. distachyon and B. stacei
plot(D.sum, S.sum)
nrow(as.matrix(D.sum))
nrow(as.matrix(S.sum))
B.Km <- cbind(D.sum, S.sum)


############ Filter by minumum call number per genome
head(B.Km)

# Assign Column 4 as Sum of Both Sub-genomes Calls
# Add Columns
Total <- (B.Km[,1])+(B.Km[,2])

# Combine The columns
B.Km <- cbind(B.Km, Total)

head(B.Km)
dim(as.matrix(B.Km))

# Plot Lowest to Highest for Total Calls
plot(B.Km[,1],B.Km[,2])
Var.Min.SD <- 10000
abline(h=Var.Min.SD, lty=2)
Var.Min.DS <- (Var.Min.SD*27)/24
abline(v=Var.Min.DS, lty=2)
B.Km <- B.Km[B.Km[,3] > Var.Min.SD,]


# Minimum counts
min(B.Km[,2])
# Check Dimensions and Format
head(B.Km)
dim(B.Km)

# Normalize Data
NrmD <- (B.Km[,1])/(B.Km[,3])
NrmS <- (B.Km[,2])/(B.Km[,3])
B.Km <- cbind(B.Km, NrmD, NrmS)
# Check Dimensions and Format
head(B.Km)
dim(B.Km)

# Set axis
B.Km.Av <- B.Km[,4:5]
plot(B.Km.Av)
dim(B.Km.Av)
head(B.Km.Av)
min(B.Km.Av)
max(B.Km.Av)
B.Km.Av <- B.Km.Av[complete.cases(B.Km.Av[,1:2]),]

# Desired k value for kmeans
k=3

# Make Kmeans Plot
cl <- kmeans(B.Km.Av, k, iter.max=1000, nstart=10)
palette(rainbow(k))
plot(B.Km[,1], B.Km[,2],col=cl$cluster, pch=16)

# Amend Matrix with Cluster IDs
Clust <- cl$cluster
B.Km <- cbind(B.Km, Clust)
min(Clust)
max(Clust)
head(B.Km)

write.table(B.Km, "Hybridum_HighConfident_and_Filtered.txt", sep="	")


# Check Dimensions and Format
head(B.Km)
dim(B.Km)

# Filter to Final Samples
B.Km <- B.Km[B.Km[,6] > 2,]

# Check Plot
plot(B.Km[,1], B.Km[,2])

rownames(B.Km)


St.check <- match(names(St), rownames(B.Km))

rownames(St) <- colnames(S)
head(S)

St <- t(St)
head(St)



ST <- cbind(St.check, St)
head(ST)
dim(ST)

# Filter Meta Data by Samples that didn't match
STF <- ST[,complete.cases(ST[,1])]
dim(STF)




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
png(file="Hybridum_Dgenome_Raw_SNP_DataImage.png",width = 3072, height = 3*768)
image(1:nrow(g),1:ncol(g),as.matrix(g),xlab = "SNPs", ylab = "samples", 
       main = paste(c("SNPs","samples"),dim(g)) )
dev.off()

################################################################################################
############################### Start Filtering ################################################
################################################################################################

# Stacei subgenome
pdf(file = "Hybridum_Sgenome_Histogram_SampleCalls.pdf",width = 30, height = 30)
######## Identify the number of genotype calls by marker
samplecalls.S<-apply(S,2, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls.S,breaks=ncol(S))
samp.thres.S <- 10000
abline(v=samp.thres.S,lty=2)
dev.off()
mean(samplecalls.S)

# Distachyon subgenome
pdf(file = "Hybridum_Dgenome_Histogram_SampleCalls.pdf",width = 30, height = 30)
######## Identify the number of genotype calls by marker
samplecalls.D<-apply(D,2, function(x) sum(!is.na(x)));
#look at distribution and try to determine a threshold for markers
hist(samplecalls.D,breaks=ncol(D))
samp.thres.D <- 11803
abline(v=samp.thres.D,lty=2)
dev.off()
mean(samplecalls.D)

pdf(file = "Hybridum_Sgenome_Histogram_SNPsPerSample.pdf",width = 30, height = 30)
################# Test a couple thresholds
table(samplecalls.S > samp.thres.S)
S.snp <- S[,samplecalls.S > samp.thres.S]
############## Look at counts across samples
snpsPerSample.S<-apply(S.snp,1,function(x) sum(!is.na(x)));
hist(snpsPerSample.S,breaks = 100)
############## Test some sample count thresholds
snpCallThresh.S <- 508
abline(v = snpCallThresh.S,lty= 2)
table(snpsPerSample.S > snpCallThresh.S)
S.minor <- S.snp[snpsPerSample.S > snpCallThresh.S,]
dev.off()


pdf(file = "Hybridum_Dgenome_Histogram_SNPsPerSample.pdf",width = 30, height = 30)
################# Test a couple thresholds
table(samplecalls.D > samp.thres.D)
D.snp <- D[,samplecalls.D > samp.thres.D]
############## Look at counts across samples
snpsPerSample.D<-apply(D.snp,1,function(x) sum(!is.na(x)));
hist(snpsPerSample.D,breaks = 100)
############## Test some sample count thresholds
snpCallThresh.D <- 508
abline(v = snpCallThresh.D,lty= 2)
table(snpsPerSample.D > snpCallThresh.D)
D.minor <- D.snp[snpsPerSample.D > snpCallThresh.D,]
dev.off()



########## Free up memory and check dimensions
rm(S)
rm(D)
rm(S.snp)
rm(D.snp)


dim(S.minor)
dim(D.minor)
colnames(S.minor)
colnames(D.minor)


############################## Rare varients Stacei
#### 1 het
rare.var.S <- rowSums(S.minor,na.rm=T) == 1
table(rare.var.S)
S.minor <- S.minor[!rare.var.S,]

##### 1 hom
rare.var.S <- (rowSums(S.minor==2,na.rm=T) == 1) & (rowSums(S.minor,na.rm=T) == 2)
table(rare.var.S)
S.minor <- S.minor[!rare.var.S,]



############################## Rare varients Stacei
rare.var.D <- rowSums(D.minor,na.rm=T) == 1
table(rare.var.D)
D.minor <- g.minor.D[!rare.var.D,]

##### 1 hom
rare.var.D <- (rowSums(D.minor==2,na.rm=T) == 1) & (rowSums(D.minor,na.rm=T) == 2)
table(rare.var.D)
D.minor <- D.minor[!rare.var.D,]







######## Duplicate markers Stacei
pdf(file = "Hybridum_Dgenome_Histogram_HetCounts.pdf",width = 30, height = 30)
het.count.S <- rowSums(S.minor==1,na.rm=T)
hist(het.count.S,breaks=ncol(S.minor))
paraThresh.S <- 5
abline(v = paraThresh.S, lty = 2)
table(het.count.S < paraThresh.S)
dev.off()


S.minor <- S.minor[het.count.S < paraThresh.S,]
########## How many total calls? ################
table(is.na(S.minor))
dim(S.minor)

######## Duplicate markers Distachyon
pdf(file = "Hybridum_Dgenome_Histogram_HetCounts.pdf",width = 30, height = 30)
het.count.D <- rowSums(D.minor==1,na.rm=T)
hist(het.count.D,breaks=ncol(D.minor))
paraThresh.D <- 5
abline(v = paraThresh.D, lty = 2)
table(het.count.D < paraThresh.D)
dev.off()


D.minor <- D.minor[het.count.D < paraThresh.D,]
########## How many total calls? ################
table(is.na(D.minor))
dim(D.minor)



######################################################################################
############################### Merge Sub-Genomes Check ##############################
######################################################################################

match(colnames(S.minor), colnames(D.minor))
match(colnames(D.minor), colnames(S.minor))

#######################################################################
########## Equate Dimensions and Merge Filtered Sub-Genomes ###########
#######################################################################

# Check Column names
colnames(S.minor)
colnames(D.minor)

# Set n.minor data frames as matrix and check headings
x = as.matrix(S.minor)
y = as.matrix(D.minor)
head(x)
head(y)

# Match up Column names of y against x, x is dominant selector
X.names <- match(colnames(x),colnames(y))
Y.names <- match(colnames(y),colnames(x))

# Get number of column positions if matching, NA if not comparable
X.names <- cbind(X.names)
Y.names <- cbind(Y.names)

# Transform name comparisons to a row, form a column
X.names <- t(X.names)
Y.names <- t(Y.names)


# Check Dimensions of Transformed n.names matrices
dim(X.names)
dim(Y.names)
dim(x)
dim(y)

# Add n.names to x and y matricies via row bind
z <- rbind(X.names, x)
w <- rbind(Y.names, y)
head(z)
head(w)

# Extract columns by first row with NA in z and w
z <- z[,complete.cases(z[1,])]
w <- w[,complete.cases(w[1,])]
dim(w)
dim(y)
dim(z)
dim(x)
head(w)

# Check for NAs in column names lists from each matrix
ZW.check <- match(colnames(z),colnames(w))
ZW.check

# Remove first row with 
SH.minor <- z[-1,]
DH.minor <- w[-1,]

# Combine by rows the filtered Stacei and Distachyon matrices
HSD.minor <- rbind(SH.minor, DH.minor)
dim(HSD.minor)
dim(SH.minor)
dim(DH.minor)
mean(HSD.minor)


# Calculate percent of SNPs corresponding to each genome
SN <- nrow(SH.minor)
DN <- nrow(DH.minor)

# Percent Stacei = TSN
TSN <- ((SN)/(DN+SN))*100
TSN

# Percent Distachyon = TDN
TDN <- ((DN)/(DN+SN))*100
TDN

# Free up memory
rm(x)
rm(y)
rm(X.names)
rm(Y.names)
rm(S.minor)
rm(D.minor)

########### minor allele frequency distribution Stacei
hist(rowSums(SH.minor,na.rm=T)/2)


########### minor allele frequency distribution Distachyon
hist(rowSums(DH.minor,na.rm=T)/2)


########### minor allele frequency distribution Hybridum
hist(rowSums(HSD.minor,na.rm=T)/2)



########### Transform Stacei Matrix
testS <- t(as.matrix(SH.minor))
testS[testS==0] <- "c"
testS[testS==1] <- NA
testS[testS==2] <- "t"

########### Transform Distachyon Matrix
testD <- t(as.matrix(DH.minor))
testD[testD==0] <- "c"
testD[testD==1] <- NA
testD[testD==2] <- "t"

########### Transform Hybridum Matrix 
testH <- t(as.matrix(HSD.minor))
testH[testH==0] <- "c"
testH[testH==1] <- NA
testH[testH==2] <- "t"




# Create Priciple Coordinate Vectors For Stacei subgenome
hap.genotypes.bin.S <- as.DNAbin(testS)
hap.gene.dist.S <- dist.dna(hap.genotypes.bin.S, model = "N" , pairwise.deletion=T)
swil.SH <- cmdscale(hap.gene.dist.S, k = 3)
swil.1_SH <- data.frame(swil.SH)
tot.var.SH <- sum(sapply(swil.1_SH, var))
per.var.SH <- round((sapply(swil.1_SH, var)/tot.var.SH)*100,2)

# Create Priciple Coordinate Vectors For Distachyon subgenome
hap.genotypes.bin.D <- as.DNAbin(testD)
hap.gene.dist.D <- dist.dna(hap.genotypes.bin.D, model = "N" , pairwise.deletion=T)
swil.DH <- cmdscale(hap.gene.dist.D, k = 3)
swil.1_DH <- data.frame(swil.DH)
tot.var.DH <- sum(sapply(swil.1_DH, var))
per.var.DH <- round((sapply(swil.1_DH, var)/tot.var.DH)*100,2)


# Create Priciple Coordinate Vectors For Hybridum subgenome
hap.genotypes.bin.H <- as.DNAbin(testH)
hap.gene.dist.H <- dist.dna(hap.genotypes.bin.H, model = "N" , pairwise.deletion=T)
swil.H <- cmdscale(hap.gene.dist.H, k = 3)
swil.1_H <- data.frame(swil.H)
tot.var.H <- sum(sapply(swil.1_H, var))
per.var.H <- round((sapply(swil.1_H, var)/tot.var.H)*100,2)
dim(swil.H)
dim(HSD.minor)

# Clear Memory of unnecessary variables
rm(hap.genotypes.bin.H)
rm(hap.genotypes.bin.S)
rm(hap.genotypes.bin.D)

################################################################################################
############################## Attached Meta Data to sub-genomes ###############################
################################################################################################

################### Read in Meta Data 
Meta.Data <- read.table("Brachy_META_NO_NAs.txt", header=T, sep="	")
rownames(Meta.Data) <- Meta.Data[,1]
dim(Meta.Data)
head(Meta.Data)

# Match Names
MD.H <- match(rownames(Meta.Data), names(H.minor), incomparables=NULL)
MD.S <- match(rownames(Meta.Data), colnames(SH.minor), incomparables=NULL)
MD.D <- match(rownames(Meta.Data), colnames(DH.minor), incomparables=NULL)
MD.SCo <- match(rownames(Meta.Data), rownames(swil.SH), incomparables=NULL)
MD.DCo <- match(rownames(Meta.Data), rownames(swil.DH), incomparables=NULL)



# Combine H, S, and D Meta Data
MD.HSD <- cbind(Meta.Data, MD.H, MD.S, MD.D)

# Combine S and D Meta Data
MD.SD <- cbind(Meta.Data, MD.S, MD.D)
MD.SD.Co <- cbind(Meta.Data, MD.SCo, MD.DCo)


# Check for names that don't match and edit if necessary
MD.namecheck <- match(names(DH.minor), rownames(Meta.Data))
MD.namecheck <- cbind(names(DH.minor), MD.namecheck)
MD.namecheck

# Check Dimensions and Names
head(MD.SD)
dim(MD.SD)
head(MD.SD.Co)
dim(MD.SD.Co)



# Filter Meta Data by Samples that didn't match
MD.SD <- MD.SD[complete.cases(MD.SD$MD.S),]
MD.SD.Co <- MD.SD.Co[complete.cases(MD.SD.Co$MD.S),]



# Remove duplicate information form memory
rm(MD.H)
rm(MD.S)
rm(MD.D)
rm(Meta.Data)

# Check the dimensions of new Meta Data Matrix
head(MD.HSD)
dim(MD.SD)


########################################################################################
############ Create Pairwise Distance for Dendrograms of each genome ###################
########################################################################################

# Hybridum
snp.1cor2.H <- snp.1cor2.dist <- as.dist(1-cor(HSD.minor,use = "pairwise.complete.obs")^2)

# Stacei
snp.1cor2.S <- snp.1cor2.dist <- as.dist(1-cor(SH.minor,use = "pairwise.complete.obs")^2)

# Distachyon
snp.1cor2.D <- snp.1cor2.dist <- as.dist(1-cor(DH.minor,use = "pairwise.complete.obs")^2)


################################################################################################
############################## Create Visuals Post Filtering ###################################
################################################################################################

######################## Plot Stacei PCoA to Distachyon PCoA
pdf(file = "Bhybridum_SxDygenome_ColourbyPlate_PCoA.pdf",width = 30, height = 30)
palette(rainbow(7))
plot(swil.DH[,1], swil.SH[,1], col=MD.SD$Continent, cex = 1.5, pch = MD.SD$Continent)
text(I(swil.DH[,1] + .0006), I(swil.SH[,1]), 
                            labels = colnames(DH.minor) ,
                            cex=0.2,col="Black")
legend(350,-500,paste( c("EU", "W. Asia",  "N. Africa", "Australia", "N. America", 
                            "S. America", "S. Africa")), cex = 3,
							pch = c(1:7), col = c(1:7))
dev.off()

head(swil.SH)
head(swil.DH)

head(snp.1cor2.H)


###################### Principle Coordinate by Plate Hybridum ##################
pdf(file = "Hybridum_Hgenome_Plate_PCoA.pdf",width = 30, height = 30)
palette(rainbow(9))
plot(swil.H[,2]~swil.H[,1], xlab = paste("PC1 var", per.var.H[1]),
                            ylab = paste("PC2 var", per.var.H[2]),
                            main = "Brachy PCoA", col=MD.SD$Plate, cex = 1, pch = 16)
text(I(swil.H[,2] + .0006)~I(swil.H[,1]), 
                            labels = colnames(DH.minor) ,
                            cex=0.5,col="Black")
legend(700,200,paste( c("EU", "W. Asia",  "N. Africa", "Australia", "N. America", 
                            "S. America", "S. Africa")), cex = 2,
							col = c(1:7), lwd = 5)
dev.off()



###################### Principle Coordinate by Plate ##################
pdf(file = "Distachyon_Continent_PCoA.pdf",width = 30, height = 30)
plot(swil.H[,2]~swil.H[,1], xlab = paste("PC1 var", H.var[1]),
                            ylab = paste("PC2 var", H.var[2]),
                            main = "Brachy PCoA", cex = 2, col = MD.Dist$Plate)
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels = names.mat ,
                            cex=0.5,col="Black")
dev.off()


######################### PDF of Kinship, Heatmap ########################

# Stacei Sub-genome Kinship Matrix
pdf(file = "Hybridum_Sgenome_Genetic_Kinship.pdf", width = 30, height = 30)
image(1:nrow(testS),1:nrow(testS),as.matrix(hap.gene.dist.S))
dev.off()

# Distachyon Sub-genome Kinship Matrix
pdf(file = "Hybridum_Dgenome_Genetic_Kinship.pdf", width = 30, height = 30)
image(1:nrow(testD),1:nrow(testD),as.matrix(hap.gene.dist.D))
dev.off()

# Hybriudm Sub-genome Kinship Matrix
pdf(file = "Hybridum_Hgenome_Genetic_Kinship.pdf", width = 30, height = 30)
image(1:nrow(testH),1:nrow(testH),as.matrix(hap.gene.dist.H))
dev.off()




# Plot hclust object as Custer Dendrogram
pdf(file = "Hybridum_Dendro_HangVariation.pdf",width = 25, height = 15)
plot(hclust(snp.1cor2.S),cex = 0.2)
dev.off()

#####################################################################################################
############################### Call Groups and Sub-geno Types ######################################
#####################################################################################################



#################### Count Stacei Groups by Unrooted Tree ###########################
SG <- hclust(dist(snp.1cor2.S))
plot(as.phylo(SG), type = "unrooted")

#  Define Number of Major Groups
S.groups = 80
S.height=2.5
S.major=7

plot((SG), cex=0.3)
rect.hclust(SG, S.groups, border = "Red")

pdf(file = "Hybridum_Sgenome_Multi-UniqueGroups_Dendro.pdf",width = 50, height = 20)
plot((SG), cex=0.3, main = "Stacei Sub-Genome")
rect.hclust(SG, h=S.height, border = "Red")
dev.off()

pdf(file = "Hybridum_Sgenome_Multi-UniqueGroups_DendroNoHang.pdf",width = 50, height = 20)
plot((SG), cex=0.3, main = "Stacei Sub-Genome", hang = -1)
rect.hclust(SG, h=S.height, border = "Red")
dev.off()



################ Count Distachyon Groups by Unrooted Tree ###########################
DG <- hclust(dist(snp.1cor2.D))
plot(as.phylo(DG), type = "unrooted")

#  Define Number of Major Groups
D.groups = 80
D.height=2.5
D.major=10

plot((DG), cex=0.3)
rect.hclust(DG, D.groups, border = "Red")

pdf(file = "Hybridum_Dgenome_Multi-UniqueGroups_Dendro.pdf",width = 50, height = 20)
plot((DG), cex=0.3, main = "Distachyon Sub-Genome")
rect.hclust(DG, h=D.height, border = "Blue")
dev.off()


pdf(file = "Hybridum_Dgenome_Multi-UniqueGroups_DendroNoHang.pdf",width = 50, height = 20)
plot((DG), cex=0.3, main = "Distachyon Sub-Genome", hang = -1)
rect.hclust(DG, h=D.height, border = "Blue")
dev.off()



################ Count Hybridum Groups by Unrooted Tree ###########################

HG <- hclust(dist(snp.1cor2.H))
plot(as.phylo(HG), type = "unrooted")

#  Define Number of Major Groups
H.groups = 80
H.height=2.2
H.major=6

plot((HG), cex=0.3)
rect.hclust(HG, H.groups, border = "Red")

pdf(file = "Hybridum_Hgenome_Multi-UniqueGroups_Dendro.pdf",width = 50, height = 20)
plot((HG), cex=0.3, main = "Hybridum Genome")
rect.hclust(HG, h=H.height, border = "green4")
dev.off()

pdf(file = "Hybridum_Hgenome_Multi-UniqueGroups_DendroNoHang.pdf",width = 50, height = 20)
plot((HG), cex=0.3, main = "Hybridum Genome", hang = -1)
rect.hclust(HG, h=H.height, border = "green4")
dev.off()





########## Principle Coordinate Analysis With Kmeans for Groups ##############################

# Stacei subgenome
pdf(file = "Hybridum_Sgenome_K-11_PCoA.pdf",width = 30, height = 30)
cl <- kmeans(swil.SH[,1:2],11)
palette(rainbow(11))
plot(swil.SH[,1], swil.SH[,2], col=cl$cluster, pch = MD.SD$Continent)
legend(550,-1100,paste("Continent",c(1:9)),
							pch = c(1:9),lty=1)
legend(550,1200,paste("K Population",c(1:11)),
							col = c(1:11),lty=1)
dev.off()

# Distachyon subgenome
pdf(file = "Hybridum_Dgenome_K-11_PCoA.pdf",width = 30, height = 30)
cl <- kmeans(swil.DH[,1:2],10)
palette(rainbow(10))
plot(swil.DH[,1], swil.DH[,2], col=cl$cluster, pch = MD.SD$Continent)
legend(550,-1100,paste("Continent",c(1:9)),
							pch = c(1:9),lty=1)
legend(550,1200,paste("K Population",c(1:11)),
							col = c(1:11),lty=1)
dev.off()

# Hybridum subgenome
pdf(file = "Hybridum_Hgenome_kmeansS7_D10_PCoA.pdf",width = 30, height = 30)
Sl <- kmeans(swil.SH[,1:2],7)
Dl <- kmeans(swil.DH[,1:2],10)
palette(rainbow(10))
plot(swil.H[,1], swil.H[,2], pch=Sl$cluster, col = Dl$cluster, cex = 2)
text(swil.H[,1], swil.H[,2], 
                            labels = colnames(DH.minor) ,
                            cex=0.5,col="Black")
legend(550,-1100,paste("Stacei K",c(1:7)), cex = 3,
							pch = c(1:7),lty=1)
legend(550,1200,paste("Distachyon K",c(1:10)), cex = 3,
							col = c(1:10),lty=1)
dev.off()





###################### Output names of Dendrogram by Groups ######################

# Write Stacei Groups
write.hclust(SG, file="Hybridum_80_Unique_Sgenome_GroupsII.txt", prefix = "BSDendro", 
							h = NULL, k = 80, append = FALSE, dec = ",")

# Distachyon Groups
write.hclust(DG, file="Hybridum_80_Unique_Dgenome_GroupsII.txt", prefix = "BDDendro", 
							h = NULL, k = 80, append = FALSE, dec = ",")

# Hybridum Groups
write.hclust(HG, file="Hybridum_80_Unique_Hgenome_GroupsII.txt", prefix = "BHDendro", 
							h = NULL, k = 80, append = FALSE, dec = ",")



############################ By Groups, Groups set by user ########################


write.hclust(SG, file="Hybridum_Groups-11_Sgenome_Groups.txt", prefix = "BSDendro", 
							h = NULL, k = 11, append = FALSE, dec = ",")

# Distachyon Groups
write.hclust(DG, file="Hybridum_Groups-11_Dgenome_Groups.txt", prefix = "BDDendro", 
							h = NULL, k = 11, append = FALSE, dec = ",")

# Hybridum Groups
write.hclust(HG, file="Hybridum_Groups-11_Hgenome_Groups.txt", prefix = "BHDendro", 
							h = NULL, k = H.groups, append = FALSE, dec = ",")


############################# Write cluster files ##################################

########## Read in cluster files per genome and sub-genomes
# Read in Hybridum
Clust.H <- read.delim(file="Hybridum_80_Unique_Hgenome_GroupsII.txt", header = T)
head(Clust.H)
row.names(Clust.H) <- Clust.H[,1]
# Read in S Genome
Clust.S <- read.delim(file="Hybridum_80_Unique_Sgenome_GroupsII.txt", header = T)
head(Clust.S)
row.names(Clust.S) <- Clust.S[,1]
# Read in D Genome
Clust.D <- read.delim(file="Hybridum_80_Unique_Dgenome_GroupsII.txt", header = T)
head(Clust.D)
row.names(Clust.D) <- Clust.D[,1]

# Match names to meta data for each sub-genome cluster ID
clust.H.match <- match(rownames(Clust.H), rownames(MD.SD))
clust.S.match <- match(rownames(Clust.S), rownames(MD.SD))
clust.D.match <- match(rownames(Clust.D), rownames(MD.SD))

# Combine match data with cluster ID data
clust.H.bind <- cbind(clust.H.match, Clust.H)
clust.S.bind <- cbind(clust.S.match, Clust.S)
clust.D.bind <- cbind(clust.D.match, Clust.D)

# Filter out samples that don't ahve meta data
clust.H.names <- clust.H.bind[complete.cases(clust.H.bind[,1]),]
clust.S.names <- clust.S.bind[complete.cases(clust.S.bind[,1]),]
clust.D.names <- clust.D.bind[complete.cases(clust.D.bind[,1]),]

# Set simpler names for colums for modified meta data
Clust.H <- clust.H.names
Clust.S <- clust.S.names
Clust.D <- clust.D.names

# Filter out un-needed columns for cluster data
Clust.H <- Clust.H[,2:3]
Clust.S <- Clust.S[,2:3]
Clust.D <- Clust.D[,2:3]

match(rownames(Clust.H), rownames(MD.SD))
match(rownames(Clust.S), rownames(MD.SD))
match(rownames(Clust.D), rownames(MD.SD))





pdf(file = "Distachyon_G.groups_PCoA.pdf",width = 30, height = 30)
plot(swil.H[,2]~swil.H[,1], xlab = paste("PC1 var", per.var.H[1]),
                            ylab = paste("PC2 var", per.var.H[2]),
                            main = "Brachy PCoA", cex = 1, pch = Clust.H$Cluster)
text(I(swil.1[,2] + .0006)~I(swil.1[,1]), 
                            labels = names.mat ,
                            cex=0.5,col="Black")
dev.off()



######################### PDF of Dendrogram no hang By Height #######################

################################################################################################
################# Create Structure Input File for Each Sub-genome ##############################
################################################################################################

########## Transform Stacei sub-genome SNP Table
Struct.S <- t(as.matrix(SH.minor))
head(Struct.S)
dim(Struct.S)
# Set absence of marker
Struct.S[is.na(Struct.S)] <- -9
# Output Structure File
write.table(Struct.S, "StructSgenome.txt", sep="	", row.names = T, col.names = T)


########### Transform Distachyon sub-genome SNP Table
# Transform SNP Table
Struct.D <- t(as.matrix(DH.minor))
head(Struct.D)
dim(Struct.D)
# Set absence of marker
Struct.D[is.na(Struct.D)] <- -9
# Output Structure File
write.table(Struct.D, "StructDgenome.txt", sep="	", row.names = T, col.names = T)


############ Transform Hybridum genome SNP Table
# Transform SNP Table
Struct.H <- t(as.matrix(HSD.minor))
head(Struct.H)
dim(Struct.H)
# Set absence of marker
Struct.H[is.na(Struct.H)] <- -9
# Output Structure File
write.table(Struct.H, "StructHgenome.txt", sep="	", row.names = T, col.names = T)


################################################################################################
###### This is the end of the script for genotyping. This script was created as one piece ######
############### and here has been chopped into smaller chunks per each chapter. ################
################################################################################################

