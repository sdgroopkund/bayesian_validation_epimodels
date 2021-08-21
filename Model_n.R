##########################################################################################################
### RUNS THE HIV TRANSMISSION MODEL GIVEN THE PARAMETERS USING C ODE SOLVER RK4 (SEE C FOLDER)
##########################################################################################################

### In this version the extra parameters are defined outside the function. 

#######################################
####  Call C version of the model  ####
#######################################

system("R CMD SHLIB C/rlib.cpp C/model.cpp C/states.cpp")  # compile the C code 
dyn.load(paste("C/rlib", .Platform$dynlib.ext, sep=""))    # load the C shared library

## function to call the C code. Takes same inputs and returns the same outptu
sim.mod.n <- function(init.pop, param, fixed.param, out.times, stepsinyear = Nostepsinyear){  ## as defined in deSolve model version
  ## init.pop: c(S1, S2, S3, S4, I1, I2, AI, NAI)
  ## param: c(notatrisk, ART, scale, circmax, intervtime, propcircatbirth, betafit)
  
  ## varying epidemiologic model parameters
  notatrisk <- param[1]
  circbirth <- param[6]
  beta <- param[7]
  
  
  ## intervention parameters
  ART <- param[2]
  scale <- param[3]
  circmax <- param[4]
  intervtime <- param[5]
  trend<-param[8]##this is new!!!!!
  
  ##fixed param
  
  epsilonart <- fixed.param[1]
  epsiloncirc <- fixed.param[2]
  mu <- fixed.param[3]
  growth <- fixed.param[4]
  prograte <- fixed.param[5]
  infdeath <- fixed.param[6]
  ARTdeath <- fixed.param[7]
  
  
  
  ## initial values
  initS1 <- init.pop[1]
  initS2 <- init.pop[2]
  initS3 <- init.pop[3]
  initS4 <- init.pop[4]#!!!! NEW EXTRA STATE
  initI1 <- init.pop[5]
  initI2 <- init.pop[6]
  initAI <- init.pop[7]
  initNAI <- init.pop[8]#!!!! NEW EXTRA STATE
  
  out <- .Call("toymodel", out.times, stepsinyear,
               growth, mu, notatrisk, circbirth,
               beta, epsiloncirc, epsilonart,
               prograte, infdeath, ARTdeath,
               ART, scale, circmax, intervtime,trend,
               initS1, initS2, initS3, initS4, initI1, initI2, initAI, initNAI)#!!!!EXTRA STATES
  
  colnames(out) <- c("time", names(init.pop), "newinfections", "addcircatbirth", "newinitationsART", "newdeathsAIDS", "behaviourchange","births")
  
  return(out)
}
