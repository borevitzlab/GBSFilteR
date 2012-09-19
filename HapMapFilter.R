#change working directory
setwd("../pel2b")
setwd("Downloads/pel2b")

# read in numberic raw file
geno <- read.table("Pel2b.hmn.txt",header=T,na.strings=".")
geno <- read.table("HapMap.hmp.txt",header=T,na.strings=".")

#strip marker info
g <- geno[,-(1:11)] 
#image raw data
jpeg(file="raw.jpg")
image(as.matrix(g),main = "19k96samp")
dev.off()

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
wells <- gsub("([[:upper:]])0([[:digit:]])", replacement="\\1\\2", names.mat[5,])
plate.num <- names.matw[3,]

#lets test this in github
## more testingd...dfsdf


#plot(read counts, snp counts) to look at correlation
#need to get raw read count table
#filter at different snps

#download .csv from gDoc
#https://docs.google.com/spreadsheet/ccc?key=0Ap_U8cpJIbwIdDEyR2JvZFY5TEtvbzhqQU5Ea3d5ckE
pel.names <- read.csv("PeliSampleData - SpecimenWell.csv")

# not all samples in the sheet were run through sequecing
index <- match(paste(names.mat[3,],row.0),
            paste(pel.names$Plate.Number,pel.names$Well.Number) )

#build names table with extra information
names.mat <- rbind(   	names.mat,
                  as.character(pel.names$Well.Number[index]),
                  as.character(pel.names$Locality[index]),
                  as.numeric(pel.names$Lat[index]),
                  as.numeric(pel.names$Long[index]),
                  as.character(pel.names$EntityCode[index]),
                  as.character(pel.names$State[index]),
                  as.character(pel.names$Complex[index]),
                  as.character(pel.names$Plant.Number[index]) )
                  
names(g) <- paste(names.mat[10,],names.mat[7,],names.mat[11,],names.mat[13,])

##move at some stage and/or clean up.
names(g.minor)[2] <- NA #"AC Rosedale NSW.1"
names(g.minor)[45] <- NA # "I Mt Stromlo ACT.1"  


entity <- names.mat[10,]
pel.complex <- names.mat[12,]

entity <- entity[samples >600]
pel.complex <- pel.complex[samples >600]
entity[c(2,45)] <- NA

locality <- as.character(pel.names$Locality[index])[samples > 600]

## reduce to counts within localities
#g.minor.sum <- matrix(ncol = length(table(locality)),nrow = nrow(g.minor))
#for (i in 1:nrow(g.minor)) g.minor.sum[i,] <- tapply(as.numeric(g.minor[i,]),locality,sum,na.rm=T)
#locality.entity <- tapply(entity,locality, function (x) x[1])

#g.test <- g.minor
#g.test[is.na(g.minor)] <- 0
#structures <- data.frame(individuals = names(g.minor), locality = locality, entity =entity)
#test <- amova(samples = g.test[,-c(2,29,45)], structures = structures[-c(2,29,45),])

install.packages("adegenet")
library("adegenet")
pel.genid <- genind(tab=g.minor/2,pop=factor(entity=="AC"),type="codom")
pel.fstat <- pairwise.fst(pel.genid,res.type="matrix")
# use colors based on Entity Code
plot.col  <- c("blue","green","red","yellow","orange","grey","brown","black","white","purple","cyan")[as.numeric(factor(names.mat[10,]))]

gps.file <- cbind(Names=names(g),color=plot.col,Lat=names.mat[8,],Long=names.mat[9,])


library(BoSSA)
gps.frame <- data.frame(nom = names(g), lon = as.numeric(names.mat[8,]),
                              lat = as.numeric(names.mat[9,]) )[samples > 600,]

snp.euc.dist <- dist(t(g.minor))

pel.complexes <- names(table(pel.complex))
pdf("mantel.pdf")
for (i in pel.complexes){
distGeo <- distGPS(gps.frame[pel.complex == i,])
snp.euc.dist.complex <- dist(t(g.minor[,pel.complex == i,]))

plot(as.dist(distGeo), snp.euc.dist.complex , main = paste("Genetic vs Geographic Distances",i), xlab = "km" )
library(ade4)
mantel.pel <- mantel.randtest(as.dist(distGeo)/1000,snp.euc.dist.complex, nrepet = 10000)
abline(coef(lm(snp.euc.dist.complex~as.dist(distGeo))))
print(mantel.pel)
}
dev.off()

#filter to samples that worked and can be changed later
gps.file <- gps.file[samples > 600,]
write.csv(file = "gps.csv", gps.file)
#go to gpsvisualizer.com and upload file to convert to google earth

# identify the number of genotype calls by parker
sites<-apply(g,1, function(x) sum(!is.na(x)));
sites<-apply(g,1, function(x) sum(x!="N"))
#look at distribution and try to determine a threshold for markers
plot(cumsum(table(sites)))
hist(sites,breaks=96)
abline(v=20,lty=2)
abline(v=70,lty=3)

#test a couple thresholds
table(sites > 70)
table(sites > 20)

#721 SNPs types on more than 70 samples
g.snp <- g[sites > 70,]

# look at counts across samples
samples<-apply(g.snp,2,function(x) sum(!is.na(x)));
samples<-apply(g.snp,2,function(x) sum(x!="N"));
hist(samples,breaks = 40)
#test some sample count thresholds
table(samples > 400)
table(samples > 600)
gg <- g.snp[,samples > 400]
gg <- g.snp[,samples > 600] #select top 80 samples
# how many total calls?
table(is.na(gg))
image(as.matrix(gg))

install.packages("ape")
library(ape)


# major vs minor allele coding
freq <- rowSums(gg,na.rm=T)
tot <- rowSums(!is.na(gg))
table(freq > tot)
g.minor <- gg
g.minor[gg == 2] <- 0
g.minor[gg == 0] <- 2
g.minor[freq < tot,] <- gg[freq < tot,]

# paralogs?
image(g.minor==1)
samples<-apply(g.minor,2,function(x) grep("[A,C,G,T]",x));
image(gg.minor!="[A,C,G,T]")

#frequency distribution across samples
hist(colSums(g.minor==1,na.rm=T),breaks=25)

#threshold for too many 'het' calls in some samples??
which(colSums(g.minor==1,na.rm=T)> 200)

# paralogs/duplicate markers 
hist(rowSums(g.minor==1,na.rm=T),breaks=20)

#threshold
paralog.snp <- rowSums(g.minor==1,na.rm=T) > 40
g.minor <- g.minor[!paralog.snp,]

# minor allele frequency distribution
minor <- rowSums(g.minor,na.rm=T)/2
# dim(g.minor)
#[1] 618  81

#count number of successful genotype
table(!is.na(g.minor))
g.nohet <- g.minor
g.nohet[g.minor == 1] <- NA
# drop SNPs without both major allele classes

# calculate allele count tables. 
# calculate genotype count tables
# capture hardy weinburg statistic [eventually test sub populations]


install.packages("ape")
library(ape)

dna.raw.dist <- dist.dna(t(g.minor),model = "RAW")
snp.bin.dist <- dist(t(g.nohet), method = "binary")
snp.euc.dist <- dist(t(g.nohet))
snp.euc.dist <- dist(t(g.minor))
## operational
boot.euc.dist <- list()
njboot <- list()
bootstraps <- matrix(nrow = 159,ncol=100)
for (i in 1:100){
boot.euc.dist[[i]] <- dist(t(g.nohet[sample(1:nrow(g.nohet),replace=T),]))
njboot[[i]] <- NJ(boot.euc.dist[[i]])
bootstraps[,i] <- mergeTrees(njeu, njboot[[i]], TRUE)
}
bootstrap <- rowSums(bootstraps =="")

#snp.cor.dist <- as.dist(1-cor(g.nohet,use = "pairwise.complete.obs"))
#plot(snp.bin.dist,snp.cor.dist)
#hcor <- hclust(snp.cor.dist)
#hbin <- hclust(snp.bin.dist)
#njcor <- NJ(snp.cor.dist)
#njbin <- NJ(snp.bin.dist)
njeu <- NJ(snp.euc.dist)
#njboot <- NJ(boot.euc.dist)

# better to polarize by ancetral state and use full freq specta
pdf(file = "allele.freq618.81.pdf")
#hist(minor,breaks= 80, main = "Peli SNPs minor allele freq")
#image(as.matrix(g.minor),main = "genotype calls")
plot(hcor,main = "hclust cor dist",cex = 0.5)
plot(hbin,main = "hclust bin dist",cex = 0.5)
plot(njcor,main = "NJ cor dist",cex = 0.5)
plot(njbin, main = "NJ bin dist",cex=0.5)

bootstraps <- mergeTrees(njeu, njboot, TRUE)

pdf("bootstrap.peli.pdf")
plot(njeu,cex=0.4)
edgelabels(bootstrap, frame="none", adj=c(0.5, 0),cex=0.4)
dev.off()
write.csv(g.minor,file = "pel2bgeno-paralog.csv")


