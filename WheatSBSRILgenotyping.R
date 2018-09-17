setwd("Downloads/")
SBS_6B_GT <- read.csv("SBS_6B_GT_first_10k.csv",as.is=T,row.names = 1)
dim(SBS_6B_GT) # has rownames

SBS_6B_GT[,1] <- ifelse(substr((SBS_6B_GT[,1]),1,1) != substr((SBS_6B_GT[,1]),3,3),NA,SBS_6B_GT[,1])
SBS_6B_GT[,2] <- ifelse(substr((SBS_6B_GT[,2]),1,1) != substr((SBS_6B_GT[,2]),3,3),NA,SBS_6B_GT[,2])

# filter SNPs and exclude missing values. 
good.snp <- SBS_6B_GT[,1] != SBS_6B_GT[,2] & !is.na( SBS_6B_GT[,1]) & !is.na( SBS_6B_GT[,2])
length((good.snp))
table(good.snp)

#Filter out bad SNPs
SBS_6Bg <- SBS_6B_GT[good.snp,]
hist(colSums(is.na(SBS_6Bg)))

#apply transposes the matrix
phased = apply(SBS_6Bg, 1, function(r) { ifelse(r==".", NA, ifelse(r==r[1], 0, 1))})
dim(phased)
phased <- as.matrix(phased)

#which(colSums(is.na(as.matrix(d))) > 10)#
#IW_S18792 IW_S18831 
#128       165 
hist(colSums(is.na(as.matrix(d))))

which(colSums(is.na(as.matrix(d))) > 10)

d = dist((as.data.frame(phased)))
dim(as.matrix(d))
hc = hclust(d)
plot(hc, cex=0.4)
hc$order
dim(phased)
o = order(rowMeans(phased, na.rm=T))
image(phased[o,o], col=c("red", "blue"))
table()
