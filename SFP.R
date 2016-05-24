#borevitz@aquilegia$ more /data/html/naturalvariation/methods/SFP.R
## ATH1 chip version5 annotation

#setwd("c:/ath1/methods") # set to your current working directory
load("ath1V5.RData")
source("readcel.R")

### set up matrix
probe.OK <- matrix(FALSE,nr=712,nc=712)
probe.OK.chrom <- matrix(NA,nr=712,nc=712)
probe.OK.position <- matrix(NA,nr=712,nc=712)
probe.OK[cbind(ath1$xpos+1, ath1$ypos+1)] <- TRUE
probe.OK.chrom[cbind(ath1$xpos+1, ath1$ypos+1)] <- ath1$chr
probe.OK.position[cbind(ath1$xpos+1, ath1$ypos+1)] <- ath1$bpstart
probe.OKmm <- matrix(FALSE,nr=712,nc=712)
probe.OKmm[cbind(ath1$xpos+1, ath1$ypos+2)] <- TRUE
probe.OKmmm <- probe.OK
probe.OKmmm[probe.OKmm] <- TRUE
probe.OK.chrom.202882 <- probe.OK.chrom[probe.OK]
probe.OK.position.202882 <- probe.OK.position[probe.OK]
image(probe.OK) # to see the valid features
chrom.order <- order(probe.OK.chrom.202882,probe.OK.position.202882)

cel.files  <- list.files(pattern=".CEL$")
chip.names <- gsub(".CEL","",cel.files)
chip.names <- gsub("DNA","",chip.names)

## make sure that Col and Ler are 1,2,3,4,5,6.
# if you need to reorder because of your file names modify
chip.names <- chip.names[c(1,2,3,4,5,6)]
## So mutant and wild type could be alphabetically 7th and 8th
print(chip.names)

# consider bg.correct from bioconductor
# bg.correct(affybatch.object,method = "rma2") # but need to setup a mask

corrected.mean <- matrix(NA,nr=nrow(ath1),nc=length(cel.files))
for (i in 1:length(cel.files)) {
   corrected.mean[,i] <- read.cel(cel.files[i], cel.image = F,
    spatial.correct = T, median = F, array.size = 712, 
    filter.size = 51, jpeg.save = T)[[1]] #only save corrrected.pm
}
library(affy)
corrected.mean <- normalize.quantiles(corrected.mean)
colnames(corrected.mean) <- chip.names
save(corrected.mean,file = "corrected.mean.RData",compress=T)

pairs(corrected.mean[sample(seq(nrow(corrected.mean)), 6000),], pch=".", lower.panel = NULL)
round(cor(corrected.mean),3)
image(cor(corrected.mean))

### from Tusher et al PNAS 2001 SAM d score

nRef <- 3 # number of replicates for the Reference Col
nEco <- 3 # number of replicates for the Ecotype
a <- (1/nRef + 1/nEco)/(nRef +nEco - 2)

# Reference Col and Ler come need to be ordered 1:6
meanRef <- rowMeans(corrected.mean[,1:nRef])
meanEco <- rowMeans(corrected.mean[,(1:nEco)+nRef])
SSRef <- rowSums((corrected.mean[,1:nRef] - meanRef)^2)
SSEco <- rowSums((corrected.mean[,(1:nEco)+nRef] - meanEco)^2)
si <- sqrt(a*(SSRef + SSEco))
s0 <- median(si)
dstat <- (meanRef - meanEco)/(si +s0) # d statistic in SAM
save(dstat, file = "dstat.RData", compress = T)

# choose how many SFPs to call based on FDR

########  FDR by permutations
choose(6,3) # number of permutations for 6 arrays, 3 per group
samp.1 <- expand.grid(0:1,0:1,0:1,0:1,0:1,0:1) # for 6 chips 20 perms
samp.1 <- samp.1[apply(samp.1,1,sum)==3,]
samp.1 <- as.matrix(samp.1)

perm.d <- matrix(nrow = length(dstat), ncol = nrow(samp.1))
for (i in 1:nrow(samp.1)){
perm.mean <- corrected.mean[,order(samp.1[i,])]
meanRef <- rowMeans(perm.mean[,1:nRef]) # since they are ordered this way in corrected mean
meanEco <- rowMeans(perm.mean[,(1:nEco)+nRef])
SSRef <- rowSums((perm.mean[,1:nRef] - meanRef)^2)
SSEco <- rowSums((perm.mean[,(1:nEco)+nRef] - meanEco)^2)
si <- sqrt(a*(SSRef + SSEco))
# s0 <- median(si) # don't recalculate s0 in perms conservative?
perm.d[,i] <- sort((meanRef - meanEco)/(si +s0))
cat(i)
gc() # garbage can clean up memory
}
rm(perm.mean)

dstat.sort <- sort(dstat)
dstat.null <- rowMeans(perm.d)

# 1 sided SAM way to threshold
deltas <- 1:8/5
res.mat <- matrix(nrow = length(deltas),ncol = 5)
colnames(res.mat) <- c("delta","ori.data","perm.data","difference","FDR")
for (j in 1:length(deltas)) {
threshold.mat <- matrix(nrow = nrow(perm.d), ncol = nrow(samp.1))
for (i in 1:nrow(samp.1)) threshold.mat[,i] <- perm.d[,i] > dstat.null + deltas[j]
sam.tab <- cbind(table(dstat.sort > dstat.null + deltas[j]), table(factor(threshold.mat, c(FALSE,TRUE)))/nrow(samp.1))
res.mat[j,] <- c(deltas[j],sam.tab[2,],sam.tab[2,1]-sam.tab[2,2],sam.tab[2,2]/sam.tab[2,1])
cat(".")
gc()
}

print(res.mat)
#     delta ori.data perm.data difference        FDR
#[1,]   0.2    20117  19074.85    1042.15 0.94819556
#[2,]   0.4    18036  10855.95    7180.05 0.60190452
#[3,]   0.6    16781   3614.30   13166.70 0.21538049
#[4,]   0.8    15898   1755.10   14142.90 0.11039753
#[5,]   1.0    15175   1166.40   14008.60 0.07686326
#[6,]   1.2    14577    911.55   13665.45 0.06253344
#[7,]   1.4    14041    795.30   13245.70 0.05664126
#[8,]   1.6    13512    729.95   12782.05 0.05402235
## you cannot get below 5% FDR with only 20 permutations 
## since 1/20 is the original data.

SFP.num <- 15000 # is more than enough < 0.10 FDR
threshold <- quantile(dstat,(nrow(ath1) - SFP.num)/nrow(ath1))
markers <- dstat > threshold

## Bulk QTL mapping on markers

bulk <- cbind(corrected.mean[markers,], ath1[markers,])
save(bulk,file="bulk.coller.RData",compress = T)
write.table(bulk,file = "LerSFP.csv", sep = ",", quote = F, row.names = F)
