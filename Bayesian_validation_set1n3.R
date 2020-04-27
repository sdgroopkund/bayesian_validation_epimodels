################# BAYESIAN VALIDATION ####################

source(file='cluster/combo_prev/bayes_valid/function_gen.R')
source(file='cluster/combo_prev/bayes_valid/Model_n.R')
source(file='cluster/combo_prev/bayes_valid/FitModel_n.R')
source(file='cluster/combo_prev/bayes_valid/GetResults_simu.R')


library(mvtnorm)
library(GPfit)
library(MCMCpack)
library(mcmc)
library(MCMCpack)

#### Initial parameters ##########

#controls the change in beta over time
trenduncertain=0 #0 is off

nloc<-6 #the number of locations
prev13<-c(0.11,0.35,0.001,0.03,0.24,0.15) #prevalence in each location in the year 2013
InitialPopSize<-10000

#Times used in model
MaxTime<-50 #the number of years the model is run for
inttimeon<-35 #the year the intervention is turned on if the intervention is implemented
fittime<-33 #the year we have prevalence data for which we fit to
Nostepsinyear<-20 # 20 steps in a year i.e. time step of 0.05 used when solving the ODEs
stepsize<-(1/Nostepsinyear) #the stepsize 0.05
Interventioninput<-c(0.1,0.05,0.15)# Intervention input at 2015, at baseline ART=0, 1-scale=0, circmax=0
sample.size=2000

bias.pct=as.character(commandArgs(trailingOnly=T)[1])
script=as.integer(commandArgs(trailingOnly=T)[2]) #to run parallel with different seeds

## controls the variance of prior on u

prior='logitgauss.prior'
update.theta=as.logical(commandArgs(trailingOnly=T)[3])


#bias.pct<-as.character(bias.pct)
bias<-gsub("[.]", "", bias.pct)
time.list=800
MCMC.samples=5000
burn.in=3000
thin.in=1
tot.samples=5000
use.X=T
int.on=T
truncated.norm=F
scale1=0.15
scale2=0.1
c_beta=100
c_invgam=15
delta.prec=0.01

################################################################
################### Fixed model values as original (Y_R) ########
#################################################################


## the values the modelers give us. True u:

beta.m<-0.6 # rate of transmisson, bet 0-1
InitialPropInfect.m<-0.005 # the initial prevalence, bet 0-1
circatbirth.m<-0.25 # c in the model, bet 0-1 used to specify initials and also proportion entering circumcised group at birth
epsilonart.m <- 0.7 # \epsilon_A, bet 0-1
epsiloncirc.m <- 0.6 # \epsilon_c, bet 0-1
mu.m <- 1.0 / 35 # \mu_1, can consider bet 0-1
growth.m <- 0.03 # \epsilon, pop growth can consider bet 0-1
prograte.m <- 0.2 # \sigma_1, bet 0-1
infdeath.m <- 0.2  # \mu_2, bet 0-1
ARTdeath.m <- 0.05 # \mu_3, bet 0-1


## note that \tau is optimized using other parameters through the fitmodel function for each location

### The intervention details are also fixed
## \alpha, the proportion of people who receive ART is fixed
## \gamma, the force of infection is fixed
## max allowed circumcision the third intervention is also fixed



ini.parms<-c(beta.m,InitialPropInfect.m,circatbirth.m)
fixed.parmsfit<-c(epsilonart.m,epsiloncirc.m,mu.m,growth.m,prograte.m,infdeath.m,ARTdeath.m)
u.start<-c(fixed.parmsfit,ini.parms)
u.m<-u.start

## default width of theta prior
delta=NULL

## other simulation parameters
Z.input<-cbind(prev13,matrix(rep(u.start,each=length(prev13)),nrow=length(prev13))) ## appending the field data X with calibration parameters u

## the field data, the model output and the bias at field data augmented with true u
dat.true<-data.generate(nloc,prev13,u.start,time.list,Interventioninput,sample.size,bias.pct)




##############################################
#### returning values on the LOGIT scale  ####
##############################################


yM_int<-drop(dat.true$yM_star_int)
yM_bas<-drop(dat.true$yM_star_bas)
yR_int<-drop(dat.true$yR_star_int)
yR_bas<-drop(dat.true$yR_star_bas)
yF_int<-drop(dat.true$yF_star_int)
yF_bas<-drop(dat.true$yF_star_bas)
b_int<-drop(dat.true$b_int)
b_bas<-drop(dat.true$b_bas)

print(b_int)
print(yM_int)
print(yR_int)
print(yF_int)

parameters<-dat.true$pars.ini


########################################
##########  Priors for theta  ###########
#########################################

## getting the prior information
prior.u=get(prior)

### new parameter set needed for the data generation
u.n<-prior.u(u.start,delta.prec,delta)

# output at field data points augmented with new u
mod.data.int<-mod.output(nloc,prev13,u.n,time.list,Interventioninput,int.on=T)
mod.data.bas<-mod.output(nloc,prev13,u.n,time.list,Interventioninput,int.on=F)


#new model outputs and bias
yM_int.n<-drop(mod.data.int$yM_star)
yM_bas.n<-drop(mod.data.bas$yM_star)
b_int.n<-yF_int-yM_int.n
b_bas.n<-yF_bas-yM_bas.n

print(yM_int.n)
print(b_int.n)

# the new augmented data
Z.new<-cbind(prev13,matrix(rep(u.n,each=length(prev13)),nrow=length(prev13)))
X.new<-as.matrix(prev13)


if (use.X==T) {
  GP_int.x<-GP_fit(X.new,b_int.n,nug_thres=10,corr=list(type='exponential',power=2))
  GP_bas.x<-GP_fit(X.new,b_bas.n,nug_thres=10,corr=list(type='exponential',power=2))
  beta_int.x<-GP_int.x$beta
  beta_bas.x<-GP_bas.x$beta
  mu_int.x<-predict.GP(GP_int.x,X.new)$Y_hat
  mu_bas.x<-predict.GP(GP_bas.x,X.new)$Y_hat
  R_int.x<-corr_matrix(X.new,GP_int.x$beta, corr=list(type='exponential',power=2))
  R_bas.x<-corr_matrix(X.new,GP_bas.x$beta, corr=list(type='exponential',power=2))
  gamma_int.x<-5*1/GP_int.x$sig2 ## setting gamma for prior specification to =5*lambdab_hat
  gamma_bas.x<-5*1/GP_bas.x$sig2 ## setting gamma for prior specification to =5*lambdab_hat
  Sig_int.x<-GP_int.x$sig2*R_int.x
  Sig_bas.x<-GP_bas.x$sig2*R_bas.x
  
  if (int.on==T){
    gammab=gamma_int.x
    betab=beta_int.x
    mub=rep(0,nloc)
    yF=yF_int
    R<-corr_matrix(prev13,betab, corr=list(type='exponential',power=2))
    R.fixed=R
  } else {
    gammab=gamma_bas.x
    betab=beta_bas.x
    mub=rep(0,nloc)
    yF=yF_bas
    R.fixed=NULL
    
  }
  
} else {
  GP_int<-GP_fit(Z.new,b_int.n,nug_thres=10,corr=list(type='exponential',power=2))
  GP_bas<-GP_fit(Z.new,b_bas.n,nug_thres=10,corr=list(type='exponential',power=2))
  beta_int<-GP_int$beta
  beta_bas<-GP_bas$beta
  mu_int<-predict.GP(GP_int,Z.new)$Y_hat
  mu_bas<-predict.GP(GP_bas,Z.new)$Y_hat
  R_int<-corr_matrix(Z.new,GP_int$beta, corr=list(type='exponential',power=2))
  R_bas<-corr_matrix(Z.new,GP_bas$beta, corr=list(type='exponential',power=2))
  gamma_int<-5*1/GP_int$sig2 ## setting gamma for prior specification to =5*lambdab_hat
  gamma_bas<-5*1/GP_bas$sig2 ## setting gamma for prior specification to =5*lambdab_hat
  Sig_int<-GP_int$sig2*R_int
  Sig_bas<-GP_bas$sig2*R_bas
  
  
  if (int.on==T){
    gammab=gamma_int
    betab=beta_int
    mub=rep(0,nloc)
    yF=yF_int
    R<-corr_matrix(prev13,betab, corr=list(type='exponential',power=2))
    R.fixed=R
  } else {
    gammab=gamma_bas
    betab=beta_bas
    mub=rep(0,nloc)
    yF=yF_bas
    R.fixed=NULL
  }
}

print(mu_int.x)



# delta<-delta.start
# delta.prec<-sqrt(1e-3) #based on a hypothetical sample size of 1000


############################################################
### If we want to have separate c_beta values          #####
### ** currently using different values for groups **  #####
### Group 1 - 1,2,6,8,10 & Group 2 - 3,4,7             #####
### and currentky Group 3 - 9                          #####
############################################################

c_beta2=1.5*c_beta
c_beta3=2*c_beta
c_beta4=1.25*c_beta
c_beta_n=rep(c_beta,length(u.m))
c_beta_n[c(3,4,7)]=c_beta2
c_beta_n[5]=c_beta4


mcmc.run<-MCMC.ana.mh(prev13=prev13,yF=yF,mod.output=mod.output,u.m=u.m,mub=mub,gammab=gammab,
                      betab=betab,time.list=time.list,sample.size=sample.size,scale1=scale1,
                      scale2=scale2,c_beta=c_beta_n,c_invgam=c_invgam,truncated.norm=truncated.norm,
                      sequential=T, delta=delta,delta.prec=delta.prec,use.X=use.X,int.on=int.on,
                      Interventioninput=Interventioninput,MCMC.samples=MCMC.samples,burn.in=burn.in,
                      thin.in=thin.in,update.theta=update.theta,prior=prior)









