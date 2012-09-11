#change working directory
setwd("/home/kevin/UniWork/BIOL3157/Assignments/3")
setwd("./pel2b")

### Import data
# Read TASSEL output which contains SNP data
geno <- read.table("Pel2b.hmn.txt",header=T,na.strings=".")
# strip sample info by discarding columns 1-11
# Gives a <num.samples> X <num.snps> data frame, 96 x 19507 in the case of pel2b
g <- geno[,-(1:11)] 

### Image raw data
# Save (to png device)
png(file="raw.png")
# Create colour map of genotype matrix
# White cells are missing data
# Green is 0 Homozygous minor allele
# Brown is 1, Het
# Blue-grey is 2 Homozygous major allele
image(as.matrix(g), main="Raw Data (pel2b)", col=rainbow(3))
dev.off() 

### Filter data
# For each row (SNP) in the genotype, count how many samples have data on this snp
snp.sample.counts <- apply(g, 1, function(x) sum(!is.na(x)))
snp.cutoff = 78

# Plot the number of samples each snp has
plot(cumsum(table(snp.sample.counts)))
hist(snp.sample.counts,breaks=96)
abline(v=snp.cutoff)

# Save histogram which describes the cutoff point used
png(file="snp_filtering_hist.png")
hist(snp.sample.counts,breaks=96)
abline(v=snp.cutoff)
dev.off()

# How many SNPs are present after filtering (= TRUE)
table(snp.sample.counts > snp.cutoff)

# Transfer useful SNPs to new data frame
g.snp <- g[snp.sample.counts > snp.cutoff,]


### Qualtity control samples
# For each column (Samples) in the filterd SNP set, count how many SNPs this sample has
sample.snp.counts<-apply(g.snp,2,function(x) sum(!is.na(x)))
sample.cutoff = 500
# Plot graphical summary of samples
hist(sample.snp.counts, breaks=40, xlim=c(0,600))
abline(v=sample.cutoff)

png(file="sample_filtering_hist.png")
hist(sample.snp.counts, breaks=40, xlim=c(0,600))
abline(v=sample.cutoff)
dev.off()

table(sample.snp.counts > sample.cutoff)
gg <- g.snp[,sample.snp.counts > sample.cutoff] #select top 80 samples

# How much data is still missing?
table(!is.na(gg))
image(as.matrix(gg), col=rainbow(3))
png(file="data_post_filtering.png")
image(as.matrix(gg), col=rainbow(3))
dev.off()


### Calculate major and minor allele encodings, and filter paralogs
freq <- rowSums(gg, na.rm=T)
tot <- rowSums(!is.na(gg))
# How many samples have major allele encoded as 2?
table(freq > tot)

## re-encode minor allele homo = 2
g.minor <- gg
# Invert encoding
g.minor[gg == 2] <- 0
g.minor[gg == 0] <- 2
# replace encoding for samples which already had minor allele homo = 2
g.minor[freq < tot,] <- gg[freq < tot,] 

# Compute minor allele frequencies
min.allele.freq <- rowSums(g.minor, na.rm=T) / (rowSums(!is.na(g.minor)) *2)
hist(min.allele.freq, breaks=20)

### Filter for paralogs
# Excessive heterozygosity across many samples indicates probable paralogs
image(g.minor==1, col=rainbow(3))

## Filter samples for excessive hets
#across samples
hets.per.sample <- colSums(g.minor==1,na.rm=T)
hist(hets.per.sample, breaks=25)

#threshold for too many 'het' calls in some samples??
which(hets.per.sample> 200)

## Filter snps for excessive heterozygosity
hets.per.snp <- rowSums(g.minor==1, na.rm=T)
hist(hets.per.snp, breaks=20)

# filter out paralogs
paralogous.snps <- hets.per.snp >30
g.final <- g.minor[!paralogous.snps,]

# minor allele frequency distribution
final.min.allele.freq <- rowSums(g.final,na.rm=T)/ (2 * rowSums(!is.na(g.final)))

### Create name matrix. Do it down here, so we can apply this to the phylogenetic grouping
# Names gets the name of each column in g.final (the name of the sample)
names <- names(g.final)
# strsplit splits names of the samples in g on "_"
split.names <- strsplit(names ,split="_")
# unlist creates a single list from the split names
names.lst <- unlist(split.names)
# Create names matrix. nrow splits list into matrix with 5 rows
names.matrix <- matrix(names.lst, nrow=5)
names <- names.matrix[5,]

### Edits name matrix
# Substitute "A01" to "A1". "\\1" is a regex backreference to the first subexpression in parenthesis.
wells <- gsub("([[:upper:]])0([[:digit:]])", replacement="\\1\\2", names.matrix[5,])
plate.num <- names.matrix[3,]

# Import detailed sample data
pel.names <- read.csv("PeliSampleData - SpecimenWell.csv")

# not all samples in pel.names were run through sequecing
# this gets the subset of pel.names which were run
index <- match(
  paste(plate.num, wells),
  paste(pel.names$Plate.Number,pel.names$Well.Number)
)

# build names table with extra information
names.matrix <- rbind(
  names.matrix,
  as.character(pel.names$Well.Number[index]),
  as.character(pel.names$Locality[index]),
  as.numeric(pel.names$Lat[index]),
  as.numeric(pel.names$Long[index]),
  as.character(pel.names$EntityCode[index]),
  as.character(pel.names$State[index])
)

# Rename samples to "<entity.code> <collection.location> <collection.state>"
names(g.final) <- paste(names.matrix[10,],names.matrix[7,],names.matrix[11,])


### Make tree, define colour groups
tree <- hclust(as.dist(1-cor(g.minor,use = "pairwise.complete.obs")))
tree.groups <- cutree(tree, k=5)

# Assigns each phylogenetic group a colour
plot.col <- rainbow(max(as.numeric(tree.groups)))[as.numeric(tree.groups)]

# better to polarize by ancetral state and use full freq specta ????

pdf(file = "filtered.pel2b.pdf")
# Minor Allele Frequency histogram of final data
hist(final.min.allele.freq,breaks= 80, main="Minor Allele Frequencies")
# Genotype call image of final data
image(as.matrix(g.final), main="Genotype Call Data", col=rainbow(3), ylab="Samples", xlab="SNPs")
# Plot tree, including coloured boxes around different phylogenetic groups
plot(tree, cex = 0.5)  #cex scales text by 0.5
rect.hclust(tree,k=10, border=rainbow(max(as.numeric(tree.groups))), cluster=tree.groups)
dev.off()

### Export final data.
write.csv(g.final,file = "pel2b_FilteredGenotypes.csv")


### Export data to gps file
# Creates a data.frame of each sample
gps.file <- cbind(
  Names=names(g.final),
  color=plot.col,
  Lat=names.matrix[8,],
  Long=names.matrix[9,]
  )

write.csv(file="pel2b_GPS.csv", gps.file)
