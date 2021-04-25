# myfile <- "/home/nico/Documents/R/Q_seine_1950_2008.csv"
# data <- read.table(myfile,header=FALSE,sep="\t",dec=".")
# y <- data[,2]-mean(data[,2])
# tps <- data[,1]/365.25+1900
# 
# mywd <- "/home/nico/Documents/Publications/essais/Publis Aurelien/indices/indices"
# setwd(mywd)
# data <- read.table("SOI.csv",header=TRUE,sep=";",dec=".")
# y <- data[,2]
# tps <- data[,1]/12
# yloessmod <- loess(y~tps,span=1)
# trend <- predict(yloessmod)
# y <- y-trend
# #y <- ts(y,start=tps[1],end=tps[length(tps)],deltat=tps[2]-tps[1],freq=1/(tps[2]-tps[1]))
# 
# # Signal synthétique intermittent
# temp <- rep(c(rep(0,400),sin(seq(0,5*pi,length.out=200)),rep(0,400)),10)
# test <- y <- y.ini <- sin(seq(0,10*pi,length.out=10000))+temp/15
# plot(test,type="l")
# test3 <- test[1:length(test)/2]
# temp2 <- rep(c(rep(0,400),rnorm(200),rep(0,400)),10)
# test2 <- y <- y.ini <- sin(seq(0,10*pi,length.out=10000))+temp2/15
# plot(test2,type="l")
# # test <- test[1:5000]
# 
# ### Mode mixing
# tt <- seq(0, 0.1, length = 2000)
# f1 <- 160; f2 <- 60
# xt <- sin(2*pi*f1*tt)*(tt <= 0.033 | tt >= 0.067)
# # ou
# xt <- sin(2*pi*f1*tt)*(tt <= 0.05)+sin(2*pi*f2*tt)*(tt >= 0.05)
# plot(xt,type="l")

Wavean <- function(y,absc=1:length(y),detrend=1,methode="cwt/tfd",fftsourcefile='C:/Users/duranl/Documents/R/useful_f/fftrec.R',filt="none",dwtfilt="s8",sm=0,wingauss=3,nscale=150,pwstretch=0.2,titlcwt=FALSE,plot.rc=FALSE,plotmode="serial",ylim=NULL,timeconv=1,print=TRUE) {
  library(wmtsa)
  # if (is.ts(y)==TRUE) {
  #     print("L'objet à traiter ne doit pas être de classe ts")
  #     y.ini <- y <- as.numeric(y)
  # }
  
  y.ini <- ts(y)
  # y <- y.ini-mean(y.ini)
  if (is.null(absc)==TRUE) absc <- 1:length(y)
  
  if (detrend==0)  {
    y <- y
  }
  
  if (detrend==1)  {
    linfit <- lm(y~time(y))
    y <- linfit$residuals
  }
  
  if (detrend==2)  {
    yloessmod <- loess(y~absc,span=1)
    trend <- predict(yloessmod)
    y <- y-trend
  }
  
  if (methode=="dwt") {
    plot.rc=TRUE
    y <- y.ini
     #nlevels <- wavMaxLevel(n.taps=2, n.sample=length(y.ini), xform="dwt")
        #nlevels <- as.integer(round(log(length(y),base=2)))
        #???original code from Nico:
        #nlevels <- as.integer(log(length(y),base=2))
        #modif by Lea on 15/04/2019 to add one extra detail for ELia:
        nlevels <- as.integer(log(length(y),base=2))+1
        
    y2 <- seq(y[length(y)],y[length(y)],length=2^ceiling(log(length(y.ini),base=2)))
    y2[1:length(y)] <- y
    y <- ts(y2)
    y.dwt <- wavDWT(y, wavelet=dwtfilt, n.levels=nlevels, keep.series=FALSE)
    mres.y <- wavMRD(y.dwt)
    mres <- matrix(data=NA,nrow=length(y.ini),ncol=nlevels+1)
    #     print(summary(mres.y))
    if (plotmode=="parallel") {
      X11()
      op <- par(mfrow=c(nlevels+1,1)) #par(mfrow=c(nlevels+1,1),mar=c(2,3.5,0.5,1))
      ax <- F
      ylim=range(y.ini)
      col="black"
    }
    else if (plotmode=="serial") {
      op <- NULL
      ax <- T
      ylim <- range(y.ini)
      col="blue"
    }
    meanperiod <- c()
    for (i in 1:(nlevels)) {
      mres[,i] <- mres.y[i][[1]][1:length(y.ini)]
      N <- length(mres[,i])
      yfreq <- frequency(ts(mres[,i]))
      freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
      ffty <- fft(ts(mres[,i]))
      PSD <- Mod(ffty[2:floor((N/2)+1)])/max(Mod(ffty[2:floor((N/2)+1)]))
      if ((sm!=0)==TRUE) {
        PSDliss <- predict(loess(PSD~freq,span=sm))
        meanperiod[i] <- 1/freq[which.max(PSDliss)]
      }
      else if ((sm==0)==TRUE) meanperiod[i] <- 1/freq[which(PSD==1)]
      #       if (plot.rc) {
      plot(absc,mres[,i],type="l",xlab="",ylab=c("D",i),axes=ax,ylim=ylim-mean(y.ini),col=col,lwd=1.5)
      if (plotmode=="serial") lines(absc,y.ini-mean(y.ini),lwd=0.5)
      #       }
    }
    mres[,nlevels+1] <- mres.y[nlevels+1][[1]][1:length(y.ini)]
    N <- length(mres[,nlevels+1])
    yfreq <- frequency(ts(mres[,nlevels+1]))
    freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
    ffty <- fft(ts(mres[,nlevels+1]))
    PSD <- Mod(ffty[2:floor((N/2)+1)])/max(Mod(ffty[2:floor((N/2)+1)]))
    if ((sm!=0)==TRUE) {
      PSDliss <- predict(loess(PSD~freq,span=sm))
      meanperiod[nlevels+1] <- 1/freq[which(PSD==1)]
    }
    else if ((sm==0)==TRUE) meanperiod[nlevels+1] <- 1/freq[which.max(PSD)]
    #     if (plot.rc) {
    plot(absc,mres[,nlevels+1],type="l",xlab="",ylab=c("S",nlevels),ylim=ylim,axes=ax,col=col,lwd=1.5)
    if (plotmode=="serial") lines(absc,y.ini,lwd=0.5)
    axis(1,xlim=range(absc),ylim=ylim)
    #     }
    par(op)
    
    if (print==TRUE) {
      print(y.dwt$scales)
      print(meanperiod*timeconv)
    }
    results <- list(y.dwt$scales,mres,meanperiod)
    return(results)
  }
  
  if (methode=="modwt") {
    plot.rc=TRUE
    y <- y.ini
   #original code from NIco for levels:
    nlevels <- wavMaxLevel(n.taps=2, n.sample=length(y.ini), xform="modwt")
    #modif by Lea on 13/08/2019 to add one extra detail for ELia:
    nlevels <- wavMaxLevel(n.taps=2, n.sample=length(y.ini), xform="modwt") +1
    #     nlevels <- as.integer(floor(log(length(y),base=2)))
    #     y2 <- seq(y[length(y)],y[length(y)],length=2^ceiling(log(length(y.ini),base=2)))
    #     y2[1:length(y)] <- y
    #     y <- ts(y2)
    y <- ts(y)
    y.dwt <- wavMODWT(y, wavelet=dwtfilt, n.levels=nlevels, keep.series=FALSE)
    mres.y <- wavMRD(y.dwt)
    mres <- matrix(data=NA,nrow=length(y.ini),ncol=nlevels+1)
    #     print(summary(mres.y))
    if (plotmode=="parallel") {
      X11()
      op <- par(mfrow=c(nlevels+1,1),mar=c(2,3.5,0.5,1))
      ax <- F
    }
    else if (plotmode=="serial") {
      op <- NULL
      ax <- T
      ylim <- range(y)
    }
    meanperiod <- c()
    for (i in 1:(nlevels)) {
      mres[,i] <- mres.y[i][[1]][1:length(y.ini)]
      N <- length(mres[,i])
      yfreq <- frequency(ts(mres[,i]))
      freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
      ffty <- fft(ts(mres[,i]))
      PSD <- Mod(ffty[2:floor((N/2)+1)])/max(Mod(ffty[2:floor((N/2)+1)]))
      meanperiod[i] <- 1/freq[which(PSD==1)]
      #       if (plot.rc==TRUE) {
      plot(absc,mres[,i],type="l",xlab="",ylab=c("D",i),axes=ax,ylim=ylim-mean(y.ini),col="blue",lwd=1.5)
      if (plotmode=="serial") lines(y.ini-mean(y.ini),lwd=0.5)
      #       }
    }
    mres[,nlevels+1] <- mres.y[nlevels+1][[1]][1:length(y.ini)]
    N <- length(mres[,nlevels+1])
    yfreq <- frequency(ts(mres[,nlevels+1]))
    freq <- seq.int(from=yfreq/N,by=yfreq/N,length.out=floor(N/2))
    ffty <- fft(ts(mres[,nlevels+1]))
    PSD <- Mod(ffty[2:floor((N/2)+1)])/max(Mod(ffty[2:floor((N/2)+1)]))
    meanperiod[nlevels+1] <- 1/freq[which(PSD==1)]
    #     if (plot.rc==TRUE) {
    plot(absc,mres[,nlevels+1],type="l",xlab="",ylab=c("S",nlevels),axes=ax,ylim=ylim,col="blue",lwd=1.5)
    if (plotmode=="serial") lines(y.ini,lwd=0.5)
    axis(1,xlim=range(absc))
    par(op)
    #     }
    #     print(meanperiod*timeconv)
    results <- list(y.dwt$scales,mres,meanperiod)
    return(results)    
  }
  if (methode=="cwt/tfd") {
    # WMTSA CWT + reconstruction FFT
   source(fftsourcefile)
    nextp2 <- 2*(2^ceiling(log(length(y),base=2)))
    abscpad <- c(time(y),rep(0,nextp2-length(y)),rep(0,nextp2))
    ypad <- c(y,rep(0,length(abscpad)-length(y)))
    y.cwt <- wavCWT(ypad,scale.range=deltat(ypad)*c(1, length(ypad)*0.67),n.scale=nscale,wavelet="morlet",shift=6)
    X11();plot(y.cwt, series=FALSE,grid.size=default,power.stretch=pwstretch,zoom=c(1,length(y),1,nextp2))
    rec <- NA
    #     if (plot.rc==TRUE) {
    Ds <- locator(2)$y[1:2]
    FPer <- 2^Ds*1.03
    title(c(round(FPer[1],2),round(FPer[2],2)),line=-2)
    abline(h=Ds[1],lty="dashed");abline(h=Ds[2],lty="dashed")
    ans <- readline("Poursuivre ou recommencer (p/r)?")
    while (ans=="r") {
      #     ans <- readline("Poursuivre ou recommencer (p/r)?")
      dev.off()
      X11();plot(y.cwt, series=FALSE,grid.size=default,power.stretch=pwstretch,zoom=c(1,length(y),1,nextp2))
      Ds <- locator(2)$y[1:2]
      FPer <- 2^Ds*1.03
      title(c(round(FPer[1],2),round(FPer[2],2)),line=-2)
      abline(h=Ds[1],lty="dashed");abline(h=Ds[2],lty="dashed")
      ans <- readline("Poursuivre ou recommencer (p/r)?")
      if (ans=="p") break
    }
    if (plot.rc==TRUE) {  
      rec <- fftrec(absc=abscpad,ypad,freqband=c(1/FPer[1],1/FPer[2]),detrend=FALSE,filt=filt,wingauss=wingauss)[1:length(y)]
      plot(y.ini,type="l");lines(rec+mean(y.ini),col="blue",lwd=2)
      if (titlcwt==TRUE) title(c(round(FPer[1],2),round(FPer[2],2)),line=-2)
    }
    results <- list(rec,round(FPer[1],2),round(FPer[2],2))
    return(results)
    #     }
    # par(ask=TRUE)
    # rc <- invisible(rec)
  }
}
