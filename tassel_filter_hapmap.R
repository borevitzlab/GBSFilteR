#change working directory
args <- commandArgs(trailingOnly=T)
# args[1] = "pel2b_crosscheck.hmnfas.csv"
# args[2] = "PeliSampleData - SpecimenWell.csv"
# args[3] = "pel2b_crosscheck"

### Set Cutoffs
#This is at the start of the file to make it easy to adjust them
snp.cutoff = 78 # SNP must be present in this many samples (78)
sample.cutoff = 500 # Sample must have data on at least this many SNPS (500)
paralog.cutoff = 35 # If snps have more than this many heterozygotes (35)

### Import data
# Read TASSEL output which contains SNP data
geno <- read.csv(args[1],header=T,na.strings=".")
# strip sample info by discarding columns 1-12
# Gives a <num.samples> X <num.snps> data frame
g <- geno[,-(1:12)] 
# Save tag seqs in a vector
g.seqtags <- as.character(geno[,2])


### Image raw data
# Save (to png device)
png(file=paste(args[3],"raw.png",sep="."))
# Create colour map of genotype matrix, White cells are missing data
image(as.matrix(g), main="Raw Data", col=rainbow(3))
dev.off() 


### Filter data
# For each row (SNP) in the genotype, count how many samples have data on this snp
snp.sample.counts <- apply(g, 1, function(x) sum(!is.na(x)))

# Plot the number of samples each snp has
# DISABLED TO ALLOW SCRIPTING
# plot(cumsum(table(snp.sample.counts)))
# hist(snp.sample.counts,breaks=96)
# abline(v=snp.cutoff)

# Save histogram which describes the cutoff point used
png(file=paste(args[3], "snp_filtering_hist.png", sep="."))
hist(snp.sample.counts,breaks=96)
abline(v=snp.cutoff)
dev.off()

# How many SNP sites are present after filtering (= TRUE)
table(snp.sample.counts > snp.cutoff)

# Transfer useful SNP sites to new data frame
g.snp <- g[snp.sample.counts > snp.cutoff,]

png(file=paste(args[3], "data_post_SNP_filt.png", sep="."))
image(as.matrix(g.snp), col=rainbow(3), main="Data After SNP filtering")
dev.off()



### Qualtity control samples
# For each column (Samples) in the filterd SNP set, count how many SNPs this sample has
sample.snp.counts<-apply(g.snp,2,function(x) sum(!is.na(x)))

# Plot graphical summary of samples
# DISABLED TO ALLOW SCRIPTING
# hist(sample.snp.counts, breaks=40, xlim=c(0,600))
# abline(v=sample.cutoff)

png(file=paste(args[3], "sample_filtering_hist.png", sep="."))
hist(sample.snp.counts, breaks=40, xlim=c(0,600))
abline(v=sample.cutoff)
dev.off()

table(sample.snp.counts > sample.cutoff)

# perform filtering 
gg <- g.snp[,sample.snp.counts > sample.cutoff] 


# How much data is still missing?
table(!is.na(gg))
#image(as.matrix(gg), col=rainbow(3))
png(file=paste(args[3], "data_post_sample_filt.png", sep="."))
image(as.matrix(gg), col=rainbow(3), main="Data after SNP and Sample Filtering")
dev.off()


### Calculate major and minor allele encodings
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
#DISABLED TO ALLOW SCRITPING ON CLUSTER
# hist(min.allele.freq, breaks=20)

# Image data after allele correction
png(file=paste(args[3], "data_post_allele_correction.png", sep="."))
image(as.matrix(g.minor), col=rainbow(3), main="Data after allele correction")
dev.off()


### Filter for paralogs
# Excessive heterozygosity across many samples indicates probable paralogs
# Inspect this visually
#DISABLED FOR SCRIPTING
# image(g.minor==1, col=rainbow(3), main="Heterozygous Data (Checking for Paralogs)")

## Filter for excessive hets across samples
hets.per.sample <- colSums(g.minor==1,na.rm=T)
#DISABLED FOR SCRIPTING
# hist(hets.per.sample, breaks=25)

#threshold for too many 'het' calls in some samples??
which(hets.per.sample> 200)

## Filter snps for excessive heterozygosity
hets.per.snp <- rowSums(g.minor==1, na.rm=T)
#DISABLED FOR SCRIPTING
# hist(hets.per.snp, breaks=20)

# filter out paralogs
paralogous.snps <- hets.per.snp > paralog.cutoff
g.final <- g.minor[!paralogous.snps,]

# minor allele frequency distribution
final.min.allele.freq <- rowSums(g.final,na.rm=T)/ (2 * rowSums(!is.na(g.final)))
#DISABLED TO ALLOW SCRIPTING
# hist(final.min.allele.freq, breaks=20)

# How much data is still missing?
table(!is.na(g.final))
# Image this
png(file=paste(args[3], "final_data.png", sep="."))
image(
  as.matrix(g.final), 
  sub=paste("sample_cutoff=", sample.cutoff, "snp_cutoff=", snp.cutoff), 
  col=rainbow(3),
  main=paste("Final Data for dataset", args[3])
)
dev.off()


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
samples.names <- read.csv(args[2])

# not all samples in samples.names were run through sequecing
# this gets the subset of samples.names which were run
index <- match(
  paste(plate.num, wells),
  paste(samples.names$Plate.Number,samples.names$Well.Number)
)

# build names table with extra information
names.matrix <- rbind(
  names.matrix,
  as.character(samples.names$Well.Number[index]),
  as.character(samples.names$Locality[index]),
  as.numeric(samples.names$Lat[index]),
  as.numeric(samples.names$Long[index]),
  as.character(samples.names$EntityCode[index]),
  as.character(samples.names$State[index])
)

# Rename samples to "<entity.code> <collection.location> <collection.state>"
names(g.final) <- paste(names.matrix[10,],names.matrix[7,],names.matrix[11,])


### Export final filtered data.
g.output <- g.final
g.output$tagseqs <- g.seqtags[as.numeric(row.names(g.output))]
write.csv(g.output,file = paste(args[3],"filtered.csv", sep="."))


# Minor Allele Frequency histogram of final data
png(file=paste(args[3], "minor.allel.freqs.png", sep="."))
hist(final.min.allele.freq,breaks=50, main="Minor Allele Frequencies")
dev.off()

# Genotype call image of final data
png(file=paste(args[3], "final_data.png", sep="."))
image(
  as.matrix(g.final), 
  main="Genotype Call Data", 
  sub=paste("sample_cutoff=", sample.cutoff, "snp_cutoff=", snp.cutoff),
  col=rainbow(3), 
  ylab="Samples", 
  xlab="SNPs"
)
dev.off()


### Print final summary
final.dim <- dim(g.final)
final.dat <- as.numeric(table(!is.na(g.final)))
print(paste("Final data matrix has", final.dim[1], "SNP sites in", final.dim[2], "samples"))
print(paste("Final data matrix has", final.dat[1], "pieces of missing data, and", final.dat[2], "pieces of data present"))

