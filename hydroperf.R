#Function to calculate the performance of hydrological models
#Klingon-Gupta Efficency (KGE-statistic) dimensionless
#with bias, error and variability
#NSE and variations of NSE (NSE of ln and sqrt)

#INPUT: 
#qobs= observed discharge in m3/s
#qsim=simulated dichqrge in m3/s
#interp = logical if need to interpolate the time series

#OUPUT:
#list
#$data =cbind(qobs,qsim),
#$beta=beta,
#$r=r,
#$gamma=gamma, 
#$KGE=KGEcoeff,
#$NSE=NSE,
#$NSEln=NSEln,
#$NSEroot=NSEroot,
#$msle=msle,
#$d=dr,
#$ve=ve


#necessary packages (path to modify)
library("zoo", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("hydroGOF", lib.loc="C:/Program Files/R/R-3.5.3/library")



hydroperf<-function(qobs,qsim,interp=F){

  qobs<-qobs
  qsim<-qsim

  
#Overall performance : KGE (Gupta et al.(2009)), NSE (Nash and Sutcliffe, 1970), Index of agreement (Wilmott, 1981)
 
  #KGE------------------------------------------------
   
  #mu (µ) is the mean runoff in m3/s 
  #(mus for simulated discharge, muo for observed discharge)
  mus<-mean(qsim)
  muo<- mean(qobs)
  
  #sigma is the standard deviation of the time series
  #sigmas, sigmao are respectiveley sd for simulated and observed discharge
  sigmas<-sd(qsim)
  sigmao<-sd(qobs)
  

#beta is the biais ratio (dimensionless)
  beta<-mus/muo
    
  #Gamma isthe variability ratio (dimensionless)
  #it is the ratio of covariances of simulated and observed dicharge
  gamma<-(sigmas/mus)/(sigmao/muo)
  
  #r is the correlation coefficient between simulated and observed runoff
  r<-cor(qobs,qsim)
  
  #KGE coefficient:
  KGEcoeff<-1-sqrt((r-1)^2+(beta-1)^2+(gamma-1)^2)
  
  #Nash-Sutcliffe efficiency----------------------------------
  NSE<-1-(sum((qobs-qsim)^2)/sum((qobs-mean(qobs))^2))
  
  
  #Index of agreement ----------------------------------------
  
  #mae= mean-absolute error; mad= mean absolute deviation
  #dr = revised index of agreement see Willmott 2015
  
  if (mae(qsim,qobs)<=2*mad(qsim,qobs)){
    dr<-1-mae(qsim,qobs)/(2*mad(qobs))}
  else if (mae(qsim,qobs)>2*mad(qobs)){
    dr<-2*mad(qobs)/mae(qsim,qobs)-1
  }
  
  
  #NSE of root mean squared discharges------------------------------
  NSEroot<-1-sum(((sqrt(qsim)-sqrt(qobs))^2))/sum((sqrt(qobs)-mean(sqrt(qobs)))^2)
  
  
  #Low flow performance--------------------------------
  #see below for NSE of logarithms and log-transformed Mean Squared Logarithmic Error (MSLE)
  
  #NSE of logarithms of discharges---------------------
  #useful for low flows
  NSEln<-1-(sum((log(qsim)-log(qobs))^2)/sum((log(qobs)-mean(log(qobs)))^2))
  
  #Log-transformed Mean Squared Logarithmic Error (MSLE)---------
  #to measure error on low flow more specifically see De Vos 2007
  
  msle<-mse(log(qsim),log(qobs))
  
  
  #Volume errors--------------------------
  
  #volume timing error (VE) see Criss et al 2008
  ve<-1-sum(abs(qsim-qobs))/sum(qobs)
  
  
  #final list gathering all indexes----------------------------
  #return a list? with qobs, qsim, all different coefficients (r, beta, gamma and KGE)
  performance<-list(data=cbind(qobs,qsim),beta=beta,r=r,gamma=gamma, KGE=KGEcoeff,NSE=NSE,NSEln=NSEln,NSEroot=NSEroot,msle=msle,d=dr,ve=ve)
  
  #remove objects not useful before ending
  rm(mus,muo,sigmas,sigmao,beta,gamma,r,KGEcoeff, NSE, NSEroot,NSEln,msle,dr,ve)
  
  
  result <- return(performance)
}