#########################################################
#### FITS THE PARAMETER TAU 
#########################################################

  fitmodel.n<- function(nloc, prev13, fixed.parmsfit=fixed.parmsfit, beta, circatbirth, InitialPropInfect){
  
  times <- seq(0,MaxTime,length=((MaxTime*Nostepsinyear)+1))  #time series for the output- 50 years, time step 0.05
  
  notatriskfit<-matrix(0,1,nloc)
  
  #Optimization routine- least squares now fitting all the locations at the same time
  #parameter vector input is the number of locations +1 :corresponding to the shared beta and then the proportion not at risk in each location
  
  lsq <- function(v){
    #the parameters we are fitting
    for(m in 1:nloc){
      notatriskfit[1,m]<-v[m]
    }
    
    #holding the results
    prevalence<-matrix(0,nloc,length(times)) 
    least.squares.prev=0
    
    for(m in 1:nloc){
      #specify input parameters
      
      
      parmsfit<-c(notatrisk=notatriskfit[1,m],ART=0,scale=1,circmax=0,inttime=(MaxTime+10),propcircbirth=circatbirth,betafit=beta,trend=trenduncertain)#intervention not turned on 
      
      
      #and the initial values !!!!! must be altered, more states
      xstart <- c(S1=(notatriskfit[1,m]*circatbirth*(InitialPopSize*(1-InitialPropInfect))),S2=((1-notatriskfit[1,m])*circatbirth*(InitialPopSize*(1-InitialPropInfect))),S3=((1-notatriskfit[1,m])*(1-circatbirth)*(InitialPopSize*(1-InitialPropInfect))),S4=(notatriskfit[1,m]*(1-circatbirth)*(InitialPopSize*(1-InitialPropInfect))),I1=(InitialPopSize*InitialPropInfect*0.5),I2=(InitialPopSize*InitialPropInfect*0.5),AI=0,NAI=0) 
      
      #run the model
      output <- sim.mod.n(xstart, parmsfit, fixed.parmsfit, times)
      
      #calculate the prevalence at each time step!!!!
      for(y in 1:length(times)){
        prevalence[1,y]<-(sum(output[y,6:9])/sum(output[y,2:9])) #!!!!!! MORE INFECTION STATES
        
        #compare modelled prevalence and input data for all time points (note only one data point at the moment)
        if(((is.nan(prevalence[1,y]))==FALSE)&&((is.na(Observedprevalence[m,y]))==FALSE)){#confirm we do have a comparison value
          least.squares.prev=least.squares.prev+((Observedprevalence[m,y]-prevalence[1,y])^2)
          #print(least.squares.prev)
        }
      }
      
    } #location loop
    
    least.squares=least.squares.prev
    return(least.squares)  
  }
  
  
  
  
  #needed for model fitting
  Noparams<-nloc #we are fitting as a set, beta and  the proportion entering the model who are not at risk in each location
  Noobserv<-1 #we are fitting to the input prevalence in year 2013 in each location
  param<-matrix(0,1,Noparams)#matrix to hold the fitted parameter values
  Observedprevalence<-matrix(NA,nloc,length(times))#table so can increase number of 'observations' in future work
  #ensure NA if no data is available for that time step.
  
  
  
  v<-c(rep(0,(nloc)))
  #initials
  for(l in 1:(nloc)){
    v[l]<-0.5 #vector of parameters to be fit- initial values
  }
  for(l in 1:nloc){
    Observedprevalence[l,((fittime*Nostepsinyear)+1)]<-prev13[l]#note subscript here has to match the number in the outputs not the time as in the model 
  }
  
  optimise<-optim(v,lsq,lower=c(rep(0.001,(nloc))),upper=c(rep(0.85,(nloc))), method="L-BFGS-B")
  param[1,]<-optimise$par[]
  converge<-optimise$convergence
  #note if look at output here (optimise) can see additional indicators- 0 indicates successful convergence
  
  
  
  
  return(param)
}
##############################################################################
