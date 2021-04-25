#Function to calculate Klingon-Gupta Efficency (KGE-statistic) dimensionless 

#INPUT: 
#qobs= observed discharge in m3/s
#qsim=simulated dichqrge in m3/s
#interp = logical (T or F) if need to interpolate the time series

#OUPUT:




hydroperf<-function(qobs,qsim,interp=F){

#necessary packages  

  qobs<-qobs
  qsim<-qsim
  
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
  
  #adding NSE?
  
  
  
  #return a list? with qobs, qsim, all different coefficients (r, beta, gamma and KGE)
  performance<-list(data=cbind(qobs,qsim),beta=beta,r=r,gamma=gamma, KGE=KGEcoeff)
  
  
  
  #remove objects not useful before ending
  rm(mus,muo,sigmas,sigmao,beta,gamma,r,KGEcoeff)
  
  
  result <- return(performance)
}