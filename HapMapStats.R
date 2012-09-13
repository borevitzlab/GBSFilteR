t<-read.table('HapMap.hmc.txt',header=TRUE,row.names=1);
c<-t[,1:32]; #discard last 5 columns
colnames(c)<-sub('^X','',colnames(c));
sites<-apply(c,1,function(x) sum(x!='0|0'));
sites.s<-sort(sites,decreasing=TRUE);
samples<-apply(c,2,function(x) sum(x!='0|0'));
samples.s<-sort(samples,decreasing=TRUE);
good.samples<-names(samples.s[samples.s>=700]); #700 as threshold drops 4 bad samples, keeping 28
f<-c[,good.samples];
sites.f<-apply(f,1,function(x) sum(x!='0|0'));
sites.fs<-sort(sites.f,decreasing=TRUE);
good.sites<-names(sites.fs[sites.fs>=7]); #7(of 28) as threshold keeps 5902 good loci
f.f<-f[good.sites,];
write.table(f.f,file='good.hmc.txt',sep="\t",quote=FALSE,col.names=NA);
g<-read.table('HapMap.hmp.txt',header=TRUE,row.names=1);
colnames(g)<-sub('^X','',colnames(g));
g.f<-g[good.sites,c(colnames(g)[1:10],good.samples)];
write.table(g.f,file='good.hmp.txt',sep="\t",quote=FALSE,col.names=NA);
library(Cairo);
Cairo(file="samples.png",type="png",width=800,height=600,pointsize=12,bg="transparent",canvas='white',units ='px',dpi='auto');
	par(mar=c(7,4,4,2)+0.1); #to leave more room for sample name labels
	barplot(samples.s,las=2,main='Loci per sample',ylab='Number of loci',xlab='Sample');
dev.off();
Cairo(file="loci.png",type="png",width=800,height=600,pointsize=12,bg="transparent",canvas='white',units ='px',dpi='auto');
	par(mar=c(7,4,4,2)+0.1);
	plot(sites.s,type='l',main='Sample coverage by loci',ylab='Samples covered by individual locus',xlab='Loci');
dev.off();
Cairo(file="loci.filt.png",type="png",width=800,height=600,pointsize=12,bg="transparent",canvas='white',units ='px',dpi='auto');
	par(mar=c(7,4,4,2)+0.1);
	plot(sites.fs[good.sites],type='l',main='Sample coverage by filtered loci',ylab='Samples covered by individual locus',xlab='Loci',ylim=c(0,max(sites.fs)));
dev.off();