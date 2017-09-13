#-------------------------------------------------
# This program is developed for the Borevitz lab
# You may use it without any warranty
# Distribution is not permitted
# Copyright - Riyan Cheng 2016
#-------------------------------------------------

#
# This program is to process raw data
# and perform GWAS
#
# Notes:
# 1) This program runs one data file at a time
# 2) Row 1: plantid; row 2: pot number
# 3) Column 1: timestamp 
# 4) Manually select days 'ud' to reliably smooth if needed
# 5) File "*-traitcapture-db-import-Plants.csv" about design
#    e.g. "BVZ0063-traitcapture-db-import-Plants.csv"
#    must immediately under the experiment folder e.g. "BVZ0063"

###################
# basic functions #
#############################################################
strf<- function(strList, pos=1){
   strTmp<- NULL
   for(i in 1:length(strList)){
      strTmp<- c(strTmp,strList[[i]][pos])
   }
   strTmp
}

# smoothing span
nf<- function(n){
   5+(30-4)/(40-5)*(n-5)
}

# remove outliers
drop1f<- function(dat, min.n=5, span=NA){
# dat: data.frame(x=, y=)
# min.n: min no. of useable observations
# span: pass to 'loess'
   drop<- NA
   predicted<- NA
   if(sum(!is.na(dat$y)) >= min.n){
      n<- sum(!is.na(dat$y))
      if(missing(span) || is.na(span))
         spn<- nf(n)/n
      lf<- loess(y ~ x, data=dat, na.action=na.exclude, degree=2, span=spn,
         control=loess.control(surface="direct"))
      rm(n)
      lfp<- predict(lf, dat$x)

      rss<- Inf
      cv<- 0
      for(i in 1:length(dat$y)){
         if(is.na(dat$y[i])) next

         .dtt<- dat; .dtt$y[i]<- NA
         .n<- sum(!is.na(.dtt$y))
         if(missing(span) || is.na(span))
            spn<- nf(.n)/.n
         .lf<- loess(y ~ x, data=.dtt, na.action=na.exclude, degree=2, span=spn,
               control=loess.control(surface="direct"))
         .lfp<- predict(.lf, .dtt$x)
         .sd<- sd(dat$y-.lfp, na.rm=TRUE) #too stringent to use .dtt
         .cv<- abs(dat$y-.lfp)[i]/(.sd + sqrt(.Machine$double.eps))
         if(FALSE){
            # large but minimize RSS
            .sd<- sd(.dtt$y-.lfp, na.rm=TRUE)
            if(.cv > qnorm(1-0.05/.n) && .sd < rss){
               lfp<- .lfp
               drop<- i
               rss<- .sd
            }
         }else{
            if(.cv > max(cv,qnorm(1-0.05/.n))){
               lfp<- .lfp
               drop<- i
               cv<- .cv
            }
         }

         rm(.n)
      }
   }
   if(!is.na(drop)){
      dat$y[drop]<- NA
      predicted<- lfp
   }
   list(data=dat, predicted=predicted, drop=drop)
}

# fit loess per day and then select the max among middle 50%
procDayDataf<- function(pdat, dl=10, min.n=0, span=NA){
# pdat: time points by row, samples by column
#    except the first column timestamp
# dl: extract day string with length of dl
# min.n: min number of observations in a day; or, no data
   itv<- table(substring(pdat$time,15,16))
      itv<- names(itv)[itv > mean(itv)/3]
      itv<- sort(as.numeric(itv),decreasing=FALSE)
      itv<- min(diff(itv))
      itv<- seq(0,60-1e-8,by=itv)
   hrs<- table(substring(pdat$time,12,13))
      hrs<- names(hrs)[hrs > mean(hrs)/3]
      hrs<- sort(as.numeric(hrs),decreasing=FALSE)
   nr<- length(hrs)
   tt<- c(paste(sprintf("%02d",min(hrs)),"_00",sep=""),paste(sprintf("%02d",max(hrs)),"_30",sep=""))
   tm<- paste(sprintf("%02d",rep(hrs,rep(length(itv),nr))), sprintf("%02d",rep(itv,nr)), sep="_")
      tm<- tm[match(tt[1],tm):match(tt[2],tm)]
   day<- sapply(as.character(pdat$time), substr, 1, dl)
   ud<- table(day)
      ud<- names(ud)[ud > mean(ud)/3]
      ud<- sort(ud, decreasing=FALSE)
   dayTime<- paste(rep(ud,rep(length(tm),length(ud))), rep(tm,length(ud)), sep="_")
   dt<- matrix(NA, nrow=length(dayTime), ncol=ncol(pdat)-1)
   rownames(dt)<- dayTime
   colnames(dt)<- colnames(pdat)[-1]
   for(j in 2:ncol(pdat)){
      if((j-1)%%10 + 1 == 1) cat(j-1) else cat(".")
      for(i in 1:length(ud)){
         #cat(".")
         dat<- pdat[day==ud[i],]
         if(nrow(dat) < min.n){
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         tmTmp<- sapply(as.character(dat$time), substr, max(dl)+2, max(dl)+6)
         dat<- dat[match(tm,tmTmp),]
         datTmp<- data.frame(x=1:nrow(dat), y = dat[, j])
         rng<- c(min(pdat[,-1],na.rm=TRUE), max(pdat[,-1],na.rm=TRUE))

         idx1<- !is.na(datTmp$y)
         #   idx1[idx1]<- datTmp$y[idx1] > 0 # remove 0's?
         if(sum(idx1) < min.n){
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         .dat<- datTmp
            .dat$y[!idx1]<- NA
         while(TRUE){
            if(sum(!is.na(.dat$y)) < min.n){
               #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
               break
            }else{
               dp1<- drop1f(.dat, min.n=min.n, span=span)
               if(is.na(dp1$drop)) break
               .dat<- dp1$data
            }
         }
         # add points that are close
         ytmp<- abs(datTmp$y-dp1$predict)
         idx<- !(ytmp > sd(.dat$y-dp1$predict, na.rm=TRUE)*qnorm(1-0.05/2) + sqrt(.Machine$double.eps))
            idx[is.na(idx)]<- FALSE
            idx<- idx & is.na(.dat$y)
         .dat$y[idx]<- datTmp$y[idx]
         rm(ytmp, idx)

         idx2<- !is.na(.dat$y)
         idx<- idx1 & idx2
         if(sum(idx) >= min.n){
            n<- sum(!is.na(.dat$y))
            # use "interpolate" to avoid unrealistic predicted values
            lf<- loess(y ~ x, data=.dat, na.action=na.exclude, degree=2, span=nf(n)/n,
               control=loess.control(surface="interpolate"))
            lfp<- predict(lf, .dat$x)
               lfp[lfp<rng[1]]<- rng[1]
               lfp[lfp>rng[2]]<- rng[2]
            dt[match(paste(ud[i], tm, sep="_"), rownames(dt)), j-1]<- lfp
            #cat("[",j-1,",",i,"]: ", dt[704:707,134], "\n",sep="")
         }else{
            #cat(paste("[", i, ",",j,"]: ", "Less than ", min.n, " usable data points\n", sep=""))
            next
         }
         #Sys.sleep(1)
      }
   }
   cat("\n")
   invisible(dt)
}

# day data & smoothing
dtf<- function(dayDat, day, ud, dl=10){
   # one data point from each day
   ddt<- matrix(NA,nrow=length(ud),ncol=ncol(dayDat))
   rownames(ddt)<- ud
   colnames(ddt)<- colnames(dayDat)
   for(i in 1:length(ud)){
      idx<- sapply(rownames(dayDat), substr, 1, dl) == ud[i]
      for(j in 1:ncol(ddt)){
         smr<- summary(dayDat[idx,colnames(ddt)[j]])
         # 'median'; better than a time point, e.g. 12:00 pm?
         ddt[i,j]<- smr["Median"]
      }
   }

   # smooth ddt
   adt<- ddt
   for(j in 1:ncol(adt)){
      dtTmp<- data.frame(x=1:nrow(ddt), y=ddt[,j])
      n<- sum(!is.na(dtTmp$y))
      if(n < 5){
         next
      }else{
         lf<- loess(y ~ x, data=dtTmp, na.action=na.exclude, degree=2, span=nf(n)/n,
              control=loess.control(surface="interpolate"))
      }
      rm(n)
      adt[,j]<- predict(lf, dtTmp$x)
   }
   adt[adt<0]<- NA

   # growth rate
   gr<- apply(adt, 2, diff)

   # relative growth rate
   rgr<- gr/adt[-nrow(adt),]

   list(ddt=ddt, adt=adt, gr=gr, rgr=rgr)
}

# quality control
qcf<- function(pdat, missing.pr=0.25, change.pr=0.25){
# missing.pr: max missing proportion
# change.pr: min change ratio
   # remove days with too many missing data
   cn<- colnames(pdat)
   pdat<- pdat[apply(is.na(pdat),1,mean)<missing.pr,]
   # remove plants with too many missing data
   ii<- apply(is.na(pdat),2,sum)
      ii<- (1:length(ii))[ii < nrow(pdat)/2]
   length(ii)
   # remove plants with little change or the control
   dmr<- apply(apply(pdat[,ii],2,range,na.rm=TRUE),2,diff)
   ii<- ii[dmr>IQR(apply(pdat[,ii],2,median,na.rm=TRUE))*change.pr]
   pdt<- pdat[,ii]
   colnames(pdt)<- cn[ii]

   pdt
}

# prepare for myScan
datProcf<- function(pdat){
   eid<- sapply(rownames(pdat),as.character)
   uid<- unique(eid)
   pdt.<- pdat[match(uid,eid),]
   nC<- rep(1,length(uid))
   for(id in uid){
      idx<- is.element(eid,id)
      if(sum(idx) < 2) next
      tmp<- pdat[idx,]
      #if(all(is.na(tmp))) cat(id, ' ')
      pdt.[match(id,uid),]<- apply(tmp,2,mean,na.rm=TRUE)
      nC[match(id,uid)]<- nrow(tmp)
   }
   pdt<- as.data.frame(pdt.)
   colnames(pdt)<- colnames(pdat)

   list(pdat=pdt, n=nC)
}

# genome scan
myScan<- function(pdat,gdat,gmap,day,cdl){
# pdat: phenotype data
# gdat: genotype data
# gmap: genetic/physical map
# day: which column of pdat
# cdl: experimental environment
####-------------------------------------------
   if(missing(cdl)){
      pdatTmp<- datProcf(pdat)
      pdt<- pdatTmp$pdat
      pdt$y<- pdt[,day]
      idx<- !is.na(pdt$y)
      pdt<- pdt[idx,]
      n<- pdatTmp$n[idx]

      rm(pdatTmp,idx)
   }else{
      uc<- unique(cdl)
      pdt<- n<- c<- NULL
      for(u in uc){
         pdatTmp<- datProcf(pdat[cdl==u,])
         pdtTmp<- pdatTmp$pdat
         pdtTmp$y<- pdtTmp[,day]
         idx<- !is.na(pdtTmp$y)
         pdt<- rbind(pdt, pdtTmp[idx,])
         n<- c(n, pdatTmp$n[idx])
         c<- c(c, rep(u,sum(idx)))

         rm(pdatTmp,pdtTmp,idx)
      }
      cdl<- c
      rm(uc,c,u)
   }
   idx<- match(rownames(pdt),rownames(gdat))
   pdt<- pdt[!is.na(idx),]
   n<- n[!is.na(idx)]
   if(!missing(cdl))
      cdl<- cdl[!is.na(idx)]
   idx<- idx[!is.na(idx)]
   gdt<- gdat[idx,]

   cat("Day ", day, " [sample size: ",nrow(pdt),"]: ", date(),"\n", sep="")

   vc<- vector("list",5)
   for(j in 1:5){
      idx<- match(rownames(gdt),rownames(gm[[j]]$AA))
      if(any(is.na(idx))) stop("something wrong...")
      v<- list(
         AA = gm[[j]]$AA[idx,idx],
         DD = NULL,
         AD = NULL,
         HH = NULL,
         MH = NULL,
         EE = diag(1/n)
      )
      if(missing(cdl)){
         vc[[j]]<- estVC(y=pdt$y,v=v)
      }else{
         vc[[j]]<- estVC(y=pdt$y,x=cdl,v=v)
      }
   }
   #### genome scan
   lrt<- est<- NULL
   for(j in 1:5){
      idx<- is.element(colnames(gdt),gmap$snp[gmap$chr==j])
      gdtTmp<- gdt[,idx]
         gdtTmp<- as.matrix(gdtTmp)
      if(missing(cdl)){
         sc<- scanOne(y=pdt$y, gdat=gdtTmp, vc=vc[[j]], intcovar=NULL, minorGenoFreq=0.05)
      }else{
         sc<- scanOne(y=pdt$y, x=cdl, gdat=gdtTmp, vc=vc[[j]], intcovar=NULL, minorGenoFreq=0.05)
      }
      lrt<- c(lrt,sc$p)
      estTmp<- matrix(unlist(sc$par),nrow=length(sc$par),byrow=TRUE)
            rownames(estTmp)<- names(sc$par)
            colnames(estTmp)<- names(sc$par[[1]])
      est<- rbind(est, estTmp)
   }

   list( lrt = lrt, est = est)
###--------------------------------
}

#########################
# process and plot data #
########################################################################################
# create folders (if unavailable) for ouput
for(f in c("graphs","images","output")){
   if(!dir.exists(f)) dir.create(f)
   rm(f)
}

#
# data file; one at a time
#
fn<- commandArgs(trailingOnly = TRUE)[1]
#
# this file is https://raw.githubusercontent.com/borevitzlab/GBSFilteR/master/BVZ0039-GC03L-C01~fullres-area.csv
#
pheno<- strsplit(fn, "\\.")[[1]][1] # best to specify the file name with data information
   pheno<- strsplit(pheno, "\\/")[[1]]
expID<- pheno[length(pheno)-4]
fld<- paste(c(pheno[1:(length(pheno)-4)],""),collapse="/")
pheno<- pheno[length(pheno)]
#
# experiment information (best to use the same file name)
#
expInfo<- paste(fld,"_data/",expID,"-traitcapture-db-import-Plants.csv",sep="")
   expInfo<- read.csv(expInfo, check.names=FALSE)
trayNum<- strf(strsplit(as.character(expInfo$Tray),"[A-Z]"))
   trayNum<- as.integer(trayNum)
#
# input data
#
pdat<- read.csv(fn, na.string=c("Na", "NaN"), header=TRUE, check.names=FALSE, skip=0)
   pdat<- pdat[apply(!is.na(pdat),1,sum,na.rm=TRUE)>ncol(pdat)*2/3,] # remove NA rows
   pdat<- pdat[apply(pdat[,-1] > 0,1,sum,na.rm=TRUE)>ncol(pdat)*2/3,] # remove 0 rows
pid<- sapply(expInfo$PlantID[match(colnames(pdat),expInfo$"Pot")],as.character)
   pid[1]<- colnames(pdat)[1]
   pid[is.na(pid)]<- paste("Pot",colnames(pdat)[is.na(pid)],sep="")
colnames(pdat)<- pid
rm(pid)

#
# process data by fitting smoothing curves
#
dl<- 10
min.n<- 7
date()
   dayDat<- procDayDataf(pdat, dl=dl, min.n=min.n, span=NA)
date()

#
# one data point / day
#
day<- sapply(as.character(pdat$time), substr, 1, 10)
ud<- sort(unique(day), decreasing=FALSE)
# ddt, adt, gr & rgr
ddt<- dtf(dayDat, day, ud)

cntNa<- apply(is.na(ddt$ddt), 2, cumsum)

#
# save results
#
# smoothed curv over days
write.csv(ddt$adt,file=paste("output/",pheno,"-smoothed.csv",sep=""))
# growth rate
write.csv(ddt$gr,file=paste("output/",pheno,"-growthRate.csv",sep=""))

# candidates for manual checking
rownames(ddt$adt)[1:min(nrow(ddt$adt),25)]
pdt<- qcf(ddt$adt[1:min(nrow(ddt$adt),25),], missing.pr=0.25, change.pr=1/3)
ex<- expInfo[is.element(expInfo$PlantID,colnames(ddt$adt)),]
exLst<- ex[!is.element(ex$PlantID,colnames(pdt)),]
rm(pdt,ex)

write.csv(exLst, file=paste(pheno,"_checkList.csv",sep=""), row.names=FALSE)

save.image(paste(pheno,".RData",sep=""))

##########################
# plotting data
#
library(jpeg)
if(file.exists("~/WinDoc/ANU/kangaroo.jpeg")) img<- readJPEG("~/WinDoc/ANU/kangaroo.jpeg") else
   img<- NA

plotf<- function(pdat,day,dayDat,ddt,cntNa,img,type=c("pdf","jpg")){
   type<- match.arg(type)
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-Cover.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   # ---------------------------------------------
   plot.new()
   text(0.5, 0.80, "NOTE")
   text(0.5, 0.75,"************************************************************************")
   text(0.5, 0.65, paste("THIS GRAPHICALLY DISPLAYS DATA AND PREDICTION\n\n",
                         "FOR THE PURPOSE OF DIAGNOSIS",sep=""))
   if(!is.na("img")) rasterImage(img,0.35,0.30,0.65,0.50)
   text(0.5, 0.55,"************************************************************************")
   text(0.5, 0.20, format(Sys.time(), "%A, %B %d, %Y"))
   text(1, 0.023, "Warning: you are running the sofware without any warranty", cex=0.95, adj=1)
   text(1, -0.01, "_________________", adj=1)
   text(1, -0.03, paste("Copyright \uA9 R. Cheng 2015",sep=""), font=3, cex=0.75, adj=1, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()

   if(type=="jpg")
      jpeg(paste("images/",pheno,"-AllScatter.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   dTmp<- sort(unique(day), decreasing=FALSE)
      udTmp<- dTmp[(1:length(dTmp))%%3==1]
   tTmp<- rownames(dayDat)
   tmIdx<- match(tTmp, sapply(as.character(pdat$time),substr,1,16))
   rng<- range(pdat[tmIdx,-1], na.rm=TRUE)
   matplot(pdat[tmIdx,-1], ylim=rng, type="p", cex=0.25, xaxt="n",
      xlab="Time", ylab="Observed", main=pheno)
   axis(1, at=match(udTmp,day[tmIdx]), tck=-0.010, labels=FALSE)
   text(match(udTmp,day[tmIdx]), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()

   # raw data + smooth curve per day for each plant
   fmt<- (ncol(pdat)-2)%/%12+1
      fmt<- (fmt-1)%/%10+1
      fmt<- paste("%0",fmt,"d",sep="")
   if(type=="pdf") par(mfrow=c(4,3), mar=c(3.0,3.0,2.0,1), mgp=c(2,1,0))
   for(j in 2:ncol(pdat)){
      if(type=="jpg"){
         jpeg(paste("images/",pheno,"-",colnames(pdat)[j],".jpg",sep=""),
            width = 720, height = 840, quality=100, res=100)
         par(mar=c(3.0,3.0,2.0,1), mgp=c(2,1,0))
         # ---------------------------------------------
      }
      plot(1:length(tTmp), rep(NA,length(tTmp)), ylim=rng, type="n", xaxt="n", cex=0.25,
         cex.axis=0.65, xlab="Time", ylab="Observed",
         main=paste("Plant ID ", colnames(pdat)[j], sep=""))
      axis(1, at=match(udTmp,day[tmIdx]), tck=-0.015, labels=FALSE)
      text(match(udTmp,day[tmIdx]), rep(min(rng)-diff(rng)*0.125,length(udTmp)),
         labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.5, xpd=TRUE)
      for(d in dTmp){
         dIdx<- grep(d,tTmp)
         tIdx<- match(tTmp[dIdx], sapply(as.character(pdat$time),substr,1,16))
         points(dIdx, pdat[tIdx,][,j], ylim=rng, type="p", cex=0.25)
         #yTmp<- dayDat[dIdx, match(colnames(pdat)[j], colnames(dayDat))]
         # the above has trouble in case duplicate plant IDs
         yTmp<- dayDat[dIdx, j-1]
            idx<- dIdx[!is.na(yTmp)]
            yTmp<- yTmp[!is.na(yTmp)]
         lines(idx, yTmp, col=2)
      }
      if(type=="jpg") dev.off()
      # ---------------------------------------------
   }

   udTmp<- ud[(1:length(ud))%%3==1]

   #1) predicted
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-Predict1.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(ddt$ddt, na.rm=TRUE)
   matplot(ddt$ddt, type="l", xaxt="n", xlab="Day", ylab="Predicted",
      main="Prediction for All Plants (Per Day)")
   axis(1, at=match(udTmp,rownames(ddt$ddt)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(ddt$ddt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()
   #2) smooth over days (one data point for each data)
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-Predict2.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(ddt$adt, na.rm=TRUE)
   matplot(ddt$adt, type="l", xaxt="n", xlab="Day", ylab="Predicted",
      main="Smoothed Prediction for All Plants")
   axis(1, at=match(udTmp,rownames(ddt$adt)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(ddt$adt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()
   #3) remove data points with too many missing data
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-Predict3.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(ddt$ddt, na.rm=TRUE)
   idx<- apply(is.na(ddt$ddt),2,sum)
   ii<- (1:length(idx))[idx < nrow(ddt$ddt)/4]
   #& remove plants with little change or the control
   dmr<- apply(apply(ddt$ddt[,ii],2,range,na.rm=T),2,diff)
   ii<- ii[dmr>median(dmr)/10]
   matplot(ddt$ddt[,ii],type="l", xaxt="n", xlab="Day", ylab="Predicted",
      main="Prediction for Plants with Missing Proportion < 25%")
   axis(1, at=match(udTmp,rownames(ddt$ddt)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(ddt$ddt)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()
   #4) accumulative number of missing data points
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-CumulNA.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(cntNa, na.rm=TRUE)
   if(rng[2]==0) rng<- c(0,1)
   matplot(cntNa, ylim=rng, type="l", xaxt="n", xlab="Day", ylab="Cumulative Number of NAs",
      main="Missing Data Pettern for All Plants (Per Day)")
   axis(1, at=match(udTmp,rownames(cntNa)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(cntNa)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   tbl<- table(cntNa[nrow(cntNa),])
   text(rep(nrow(cntNa)-0.5,length(tbl)), as.numeric(names(tbl)), tbl, cex=0.75, col=2, pos=4, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()
   #5) growth rate
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-GrowthRate.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(ddt$gr, na.rm=TRUE)
   matplot(ddt$gr, type="l", xaxt="n", xlab="Day", ylab="Growth Rate", main="")
   axis(1, at=match(udTmp,rownames(ddt$gr)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(ddt$gr)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   if(type=="jpg") dev.off()
   #6) relative growth rate
   if(type=="jpg")
      jpeg(paste("images/",pheno,"-GrowthRateRel.jpg",sep=""), width = 720, height = 840, quality=100, res=100)
   par(mfrow=c(1,1), mar=c(3.75,3.75,3.25,1), mgp=c(2.05,1,0))
   # ---------------------------------------------
   rng<- range(ddt$rgr, na.rm=TRUE)
   matplot(ddt$rgr, type="l", xaxt="n", xlab="Day", ylab="Relative Growth Rate", main="")
   axis(1, at=match(udTmp,rownames(ddt$rgr)), tck=-0.010, labels=FALSE)
   text(match(udTmp,rownames(ddt$rgr)), rep(min(rng)-diff(rng)*0.07,length(udTmp)),
      labels=gsub("_","/",substring(udTmp,6,10)), srt=35, adj=0.75, cex=0.75, xpd=TRUE)
   # ---------------------------------------------
   if(type=="jpg") dev.off()
}

plotf(pdat=pdat, day=day, dayDat=dayDat, ddt=ddt, cntNa=cntNa, img=img, type="jpg")
if(FALSE){
   cvt<- paste("convert -compress Zip -quality 100", " images/", pheno, "_Cover.jpg", sep="")
      cvt<- paste(cvt, " images/",  pheno, "_AllScatter*.jpg", sep="")
      cvt<- paste(cvt, " images/",  pheno, "_Idv*.jpg", sep="")
      cvt<- paste(cvt, " images/",  pheno, "_Prediction_*.jpg", sep="")
      cvt<- paste(cvt, " images/",  pheno, "_CumulNA.jpg", sep="")
      cvt<- paste(cvt, " images/",  pheno, "_GrowthRate*.jpg", sep="")
      cvt<- paste(cvt, " ", pheno, "Tmp.pdf", sep="")
   system(cvt)
}

pdf(paste("graphs/",pheno,".pdf",sep=""), title=pheno, height=9)
   plotf(pdat=pdat, day=day, dayDat=dayDat, ddt=ddt, cntNa=cntNa, img=img, type="pdf")
dev.off()

q("no")

########
# GWAS #
########
library(QTLRel)

load("../../allData/geno473x250k.RData")
gmap<- phyMap; colnames(gmap)<- c("snp","chr","dist")

pdat<- qcf(ddt$adt, missing.pr=0.25, change.pr=0.25)
   pdat<- t(pdat)
rownames(pdat)<- expInfo$EcotypeID[match(rownames(pdat),expInfo$PlantID)]
days<- colnames(pdat)
lrt<- matrix(NA, nrow=nrow(gmap), ncol=ncol(pdat))
   colnames(lrt)<- days
   lrt<- cbind(gmap, lrt)
for(j in 1:ncol(pdat)){
   sc<- myScan(pdat,gdat,gmap,j)
   lrt[match(names(sc$lrt),lrt$snp),days[j]]<- sc$lrt
}
rm(j,sc,pdat)
write.csv(lrt,file=paste("output/",pheno,"_LRT.csv",sep=""),row.names=FALSE)

q("no")

#################################################
# the end #
###########

