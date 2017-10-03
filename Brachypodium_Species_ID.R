################################################################################################
## <---------------------------- Window Width for Script -----------------------------------> ##

########################################## Citation ############################################
#  Species/Ploidy ID
#  Version H.2.0
#  Authors: Jared Streich
#
########################################### CONTENTS ###########################################

########### Clear all remaining variables in R, 
########### CAREFULL... you don't ruin another script running!!! 
rm(list = ls())

################################################################################################
############################## Change Working Directory ########################################
################################################################################################
setwd("your/file/path/")

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

S.sub <- read.table("myGBSGenos_S.txt",header=T,na.strings=".")
D.sub <- read.table("myGBSGenos_D.txt",header=T,na.strings=".")

################ Check Dimensions and Headers
dim(S.sub)
dim(D.sub)

head(S.sub)
head(D.sub)

# Set SNP Names / RowNames
rownames(S.sub) <- S.sub[,1]
rownames(D.sub) <- D.sub[,1]

# Species/Sub-genomic Data
S <- S.sub[,-(1:11)] 
D <- D.sub[,-(1:11)] 
colnames(S)

################### Check Dimensions and remove extra matricies 
dim(S)
dim(D)
rm(S.sub)
rm(D.sub)


# Number of Variants in B. stacei and B. distachyon Sub-genomes
Svar <- nrow(S)
Dvar <- nrow(D)



# Format Names For Data Sheet Matching
names(S)[1] <- "ref-genome_12_12_12_12"
names.mat.S <- matrix(unlist(strsplit(names(S),split="_")),nrow=5)
names(S) <- names.mat.S[1,]

# D Genome Names Format
names(D)[1] <- "ref-genome_12_12_12_12"
names.mat.D <- matrix(unlist(strsplit(names(D),split="_")),nrow=5)
names(D) <- names.mat.D[1,]


# Match Names of S and D Genomic Data
match(names(S), names(D), nomatch = NA, incomparables = NULL)

##################################################################################################
######################### Convert calls to Presence-Absence of Markers ###########################
##################################################################################################

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
colors <- c("Blue", "Red2", "Purple")
# Plot Lowest to Highest for Total Calls
pdf(file = "Hybridum_SvsD_Subgenomes_Plot.pdf",width = 30, height = 30)
plot(B.Km[,1],B.Km[,2], pch = 16, cex = 1.4, col = colors)
Var.Min.SD <- 10000
abline(h=Var.Min.SD, lty=2)
Var.Min.DS <- (Var.Min.SD*57)/66
abline(v=Var.Min.DS, lty=2)
dev.off()


pdf(file = "Hybridum_SvsD_Subgenomes_Plot.pdf",width = 30, height = 30)
plot(B.Km[,1],B.Km[,2], pch = 16, cex = 1.4)
Var.Min.SD <- 10000
abline(h=Var.Min.SD, lty=2)
Var.Min.DS <- (Var.Min.SD*57)/66
abline(v=Var.Min.DS, lty=2)
B.Km <- B.Km[B.Km[,3] > Var.Min.SD,]
dev.off()
dim(B.Km)
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


library(e1071)
svm(B.Km.Av)
# Make Kmeans Plot
cl <- kmeans(B.Km.Av, k, iter.max=1000, nstart=10)
palette(c("Red", "Blue2", "Purple"))
plot(B.Km[,1], B.Km[,2],col=cl$cluster, pch=16)
plot(B.Km.Av,col=cl$cluster, pch=16)

clst <- cbind(cl$cluster, cl$cluster)

match(rownames(clst),rownames(B.Km))
B.filt <- cbind(clst[,1], B.Km)


# 58.8:41.1 Ratio of markers in reference, 
58-12.5
58+12.5

# Calling thresholds were based on 4:1 ratio of markers, plus bias of 8.8%
# Bias as acknowledged by divinding bias in half and adding or subtracting
# to the desired thresholds:
# Distachyon 75% or greater becomes 79%
# Stacei 25% or greater becomes 29%
# Hybridum is calcualted as 37.5 to 62.5 plus or 4%, 41.5 and 66.5


dim(B.Km.Av)
# Call Distachyon samples 79% or greater 
B.d <- B.filt
B.d <- B.d[B.d[,5] >= 0.794,]
dim(B.d)

# Call Stacei samples 29% or lower
B.s <- B.filt
B.s <- B.s[B.s[,5] <= 0.294,]
dim(B.s)

# Call Hybridum samples between 45.5% to 70.5%
B.h <- B.filt
B.h <- B.h[B.h[,5] <= 0.625,]
B.h <- B.h[B.h[,5] >= 0.415,]
dim(B.h)

# Write samples to directory
dist <- B.d
stac <- B.s
hybr <- B.h
write.table(dist, file = "dist.txt", sep = "	")
write.table(stac, file = "stac.txt", sep = "	")
write.table(hybr, file = "hybr.txt", sep = "	")


# Plot of Final Sample Calls post filtering
B.all <- rbind(B.d, B.s, B.h)
dim(B.all)


B.distachyon_markers<- B.all[,2]
B.stacei_markers <- B.all[,3]

png(file="Species_ID_StrPlt.png",width = 1000, height = 1000)
plot(B.distachyon_markers, B.stacei_markers, col = B.all[,1], pch = 16, cex = 1.8)
dev.off()


# Amend Matrix with Cluster IDs
Clust <- cl$cluster
B.Km <- cbind(B.Km, Clust)
min(Clust)
max(Clust)
head(B.Km)

write.table(B.Km, "Hybridum_HighConfident_and_Filtered.txt", sep="	")

