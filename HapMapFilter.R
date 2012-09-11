#change working directory
setwd("../pel2b")
setwd("Downloads/pel2b")

# read in numberic raw file
geno <- read.table("Pel2b.hmn.txt",header=T,na.strings=".")

#strip marker info
g <- geno[,-(1:11)] 
#image raw data
jpeg(file="raw.jpg")
image(as.matrix(g),main = "19k96samp")
dev.off()

#format names for data sheet matching
names.mat <- matrix(unlist(strsplit(names(g),split="_")),nrow=5)
row.0 <- sub("0","",names.mat[5,])
tens <- grep("[[:upper:]]1",names.mat[5,])
row.0[tens] <- names.mat[5,tens] # add back F10 

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
                  as.character(pel.names$State[index]) )
names(g) <- paste(names.mat[10,],names.mat[7,],names.mat[11,])

# use colors based on Entity Code
plot.col  <- c("blue","green","red","yellow","orange","grey","brown","black","white","purple","cyan")[as.numeric(factor(names.mat[10,]))]

gps.file <- cbind(Names=names(g),color=plot.col,Lat=names.mat[8,],Long=names.mat[9,])

#filter to samples that worked and can be changed later
gps.file <- gps.file[samples > 600,]
write.csv(file = "gps.csv", gps.file)
#go to gpsvisualizer.com and upload file to convert to google earth

# identify the number of genotype calls by parker
sites<-apply(g,1, function(x) sum(!is.na(x)));
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
hist(samples,breaks = 40)
#test some sample count thresholds
table(samples > 500)
table(samples > 600)
gg <- g.snp[,samples > 500]
gg <- g.snp[,samples > 600] #select top 80 samples
# how many total calls?
table(is.na(gg))
image(as.matrix(gg))

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

#frequency distribution across samples
hist(colSums(g.minor==1,na.rm=T),breaks=25)

#threshold for too many 'het' calls in some samples??
which(colSums(g.minor==1,na.rm=T)> 200)

# paralogs/duplicate markers 
hist(rowSums(g.minor==1,na.rm=T),breaks=20)

#threshold
paralog.snp <- rowSums(g.minor==1,na.rm=T) >60
g.minor <- g.minor[!paralog.snp,]

# minor allele frequency distribution
minor <- rowSums(g.minor,na.rm=T)/2

# better to polarize by ancetral state and use full freq specta
pdf(file = "allele.freq.pdf")
hist(minor,breaks= 80, main = "Peli SNPs minor allele freq")
image(as.matrix(g.minor),main = "genotype calls")
plot(hclust(as.dist(1-cor(g.minor,use = "pairwise.complete.obs"))),cex = 0.5)
dev.off()
write.csv(g.minor,file = "pel2bgeno-paralog.csv")
