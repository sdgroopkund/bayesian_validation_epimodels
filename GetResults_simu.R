#######################################################################
#### OUTPUT RESULTS FROM THE MODEL USING THE FITTED PARAMETERS
#######################################################################

getresults.simu<-function(param,Interventioninput,InitialPopSizeadjusted,fixed.parmsfit, beta, circatbirth, InitialPropInfect){
  
  times <- seq(0,MaxTime,length=((MaxTime*Nostepsinyear)+1))  #time series for the output- 50 years, time step 0.05
  
  QALYs<-array(0,dim=c(nloc,2,MaxTime))
  costs<-array(0,dim=c(nloc,8,MaxTime))
  
  finalresultstable<-array(0,dim=c(nloc,2,MaxTime))
  
  int<-matrix(0,3,2)
  int[1,1]=Interventioninput[1]  #on rate of ART in population
  int[1,2]=0    #off ART
  int[2,1]=(1-Interventioninput[2])  #on scaling of force of infection i.e. behaviour change
  int[2,2]=1    #off behaviour change
  int[3,1]=circatbirth+Interventioninput[3]  #on, max prop circ in pop
  int[3,2]=circatbirth #off circ 
  
  baseline<-array(0,c(length(times),15,nloc))
  interventions<-array(0,c(length(times),15,nloc))
  for(l in 1:nloc){
  
    InitialPopSize=InitialPopSizeadjusted[l]##note use scaling
    xstart <- c(S1=(param[1,l]*circatbirth*(InitialPopSize*(1-InitialPropInfect))),S2=((1-param[1,l])*circatbirth*(InitialPopSize*(1-InitialPropInfect))),S3=((1-param[1,l])*(1-circatbirth)*(InitialPopSize*(1-InitialPropInfect))),S4=(param[1,l]*(1-circatbirth)*(InitialPopSize*(1-InitialPropInfect))),I1=(InitialPopSize*InitialPropInfect*0.5),I2=(InitialPopSize*InitialPropInfect*0.5),AI=0,NAI=0) 
    
    #run the model without interventions
    parmsnoint<-c(notatrisk=param[1,l],ART=int[1,2],scale=int[2,2],circmax=int[3,2],inttime=(MaxTime+10),propcircbirth=circatbirth,betafit=beta,trend=trenduncertain)#all 'off' values and interevntion time is after the end of the model
    baseline[,,l]<- sim.mod.n(xstart, parmsnoint, fixed.parmsfit, times)  # use C version
    
    #run the model with interventions
    parmswithint<-c(notatrisk=param[1,l],ART=int[1,1],scale=int[2,1],circmax=int[3,1],inttime=inttimeon,propcircbirth=circatbirth,betafit=beta,trend=trenduncertain)#all 'on' values and intervention time is inttime
    interventions[,,l] <- sim.mod.n(xstart, parmswithint,fixed.parmsfit, times)  # use C version
    
    # Find the QALY weights
    # http://image.thelancet.com/journals/lancet/article/PIIS2214109X13701724/table?tableid=tbl3&tableidtype=table_id&sectionType=red
    # init.pop: c(S1, S2, S3, S4, I1, I2, AI, NAI)
    # Those who are susceptible have a weight of 1
    # Those who are infected currently have a weight of (1-DALY) for the relevant state
    # annual output, find the final value in that year (i.e. step 21 corresponds to the end of year 1)
    
    for(z in 1:MaxTime){
      y=(z*(1/stepsize))+1 
      QALYs[l,1,z]=(baseline[y,2,l]+baseline[y,3,l]+baseline[y,4,l]+baseline[y,5,l]+((1-0.053)*baseline[y,6,l])+((1-0.547)*baseline[y,7,l])+((1-0.053)*baseline[y,8,l])+((1-0.547)*baseline[y,9,l]))
      QALYs[l,2,z]=(interventions[y,2,l]+interventions[y,3,l]+interventions[y,4,l]+interventions[y,5,l]+((1-0.053)*interventions[y,6,l])+((1-0.547)*interventions[y,7,l])+((1-0.053)*interventions[y,8,l])+((1-0.547)*interventions[y,9,l]))
    }
    
    #First ART costs, get the number on ART at the end of each year and multiply by the annual ART cost
    for(z in 1:MaxTime){
      y=(z*(1/stepsize))+1                    #annual output, time step at end of year e.g. 21
      costs[l,1,z]=(515*baseline[y,8,l])        #art (515 is an annual cost)
      costs[l,2,z]=(515*interventions[y,8,l])   
    }
    
    #Second the circumcision costs, get the difference in the number of individuals entering the model who are circumcised (vs the number entering who are circumcised with no intervention)
    for(z in 1:MaxTime){
      for(y in 1:(1/stepsize)){
        k=((z-1)*(1/stepsize))+y+1        #this has to corresond to year definition above
        #print(z)
        #print(k)                                   #need to sum all time steps in the relevant year
        costs[l,3,z]=costs[l,3,z]+(60*baseline[k,11,l])      #need all of individuals who are circumcised due to intervention entering the model
        costs[l,4,z]=costs[l,4,z]+(60*interventions[k,11,l])
      }
    }
    
    #Third the behaviour change costs, from the total population at one point in year and assumed cost of 10 dollars per year (N.B. better represented as an annual or one-off cost?) (note this is output specifically from model, due to reasons outlined below)
    for(z in 1:MaxTime){
      k=((z-1)*(1/stepsize))+1        #this summing all values across time steps in that year                   
      costs[l,5,z]=(10*(baseline[k,14,l]))     #As this is applied to whole population, could count total pop from output but may cause problems if no indication of if the intervention is on or not
      costs[l,6,z]=(10*(interventions[k,14,l]))      #Created extra variable output from model behaviourchangecount which outputs the total pop size if behaviour change on
    }
    
    #Fourth the testing costs, based on the number of ART initiations at each time step and working back to estimate the number of tests performed in population
    for(z in 1:MaxTime){
      for(y in 1:(1/stepsize)){
        k=((z-1)*(1/stepsize))+y+1                            #need to sum all time steps relevant (i.e. all time steps within the year are summed)
        costs[l,7,z]=costs[l,7,z]+((20*baseline[k,12,l])*(sum(baseline[k,2:9,l])/sum(baseline[k,6:9,l])))      #testing (need all performed, not just annual). Costs of tests per ART initiation
        costs[l,8,z]=costs[l,8,z]+((20*interventions[k,12,l])*(sum(interventions[k,2:9,l])/sum(interventions[k,6:9,l]))) #scale up the cost as need to back calculate number tested to get those initiations
      }
    }
    
    #intervention-interventions
    for(z in 1:MaxTime){
      finalresultstable[l,1,z]<-(QALYs[l,2,z]-QALYs[l,1,z]) #intervention-baseline
      finalresultstable[l,2,z]<-(costs[l,2,z]-costs[l,1,z])+(costs[l,4,z]-costs[l,3,z])+(costs[l,6,z]-costs[l,5,z])+(costs[l,8,z]-costs[l,7,z])#intervention-baseline
      #note however that interventions costs will consistently be 0 as no intervention is applied under interventions
    }
  } 
  
  return(list(finalresultstable=finalresultstable,interventions=interventions,baseline=baseline))
}


##########################################################################
