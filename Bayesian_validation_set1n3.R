###################################################################################
### CODES TO RUN THE BAYESIAN VALIDATION FRAMEWORK (SETTINGS 1 AND 3)
###################################################################################

source(file = 'cluster/combo_prev/bayes_valid/function_gen.R')
source(file = 'cluster/combo_prev/bayes_valid/Model_n.R')
source(file = 'cluster/combo_prev/bayes_valid/FitModel_n.R')
source(file = 'cluster/combo_prev/bayes_valid/GetResults_simu.R')

library(mvtnorm)
library(GPfit)
library(MCMCpack)
library(mcmc)
library(MCMCpack)


#### 1. DECLARE INITIAL PARAMETERS

trenduncertain = 0 #Controls the change in beta over time; 0 is off (so beta is fixed)
nloc = 6 #The number of communities
prev13 = c(0.11, 0.35, 0.001, 0.03, 0.24, 0.15) #HIV prevalence at each location in the year 2013
InitialPopSize = 10000 #Initial population size at each community (1980)
MaxTime = 50 #The number of years the model is run for
inttimeon = 35 #The year the intervention is turned on if the intervention is implemented (2015)
fittime = 33 #The year we have prevalence data for which we fit to (2013)
Nostepsinyear = 20 #20 steps in a year i.e. time step of 0.05 used when solving the ODEs
stepsize = (1 / Nostepsinyear) #The step size is 0.05
Interventioninput = c(0.1, 0.05, 0.15) #Intervention input at 2015; at baseline ART=0, 1-scale=0, circmax=0
sample.size = 2000 #Total sample size for the field data
prior ='logitgauss.prior' #prior to use
time.list = 800 #Time point of interest 800 / 20 = 40 years, or in other word at the begininng of 1980 + 40 = 2020

## Simulation setting parameters
bias.pct = as.character(commandArgs(trailingOnly = T)[1]) #Bias pct introduced in model outputs to create the faulty structure scenario (parallel script argument 1)
script = as.integer(commandArgs(trailingOnly = T)[2]) #to run parallel with different seeds (parallel script argument 2)
update.theta = as.logical(commandArgs(trailingOnly = T)[3]) #whether to update priors on parameter  (parallel script argument 3)
bias = gsub("[.]", "", bias.pct)

## MCMC parameters
MCMC.samples = 5000 #Total MCMC samples
burn.in = 3000 #Total burn in
thin.in = 1 #No thinning
delta.prec = 0.01 #precision parameter for controlling the variance of the parameters logitgauss or uniform prior
delta = NULL # default width of logitgauss or uniform prior for parameters

## Metropolis Hastings hyperparameters
scale1 = 0.15 #variance of Normal distribution proposal for discrepancy at each community 
scale2 = 0.1 #variance of truncated Normal distribution proposal for gamma_b (hyperparameter for lambda_b)
c_beta = 100 #controls the parameters of the beta distribution proposal for the parameters
c_invgam = 15 #controls the parameters of the inverse gamma distribution proposal for gamma_b (hyperparameter for lambda_b)

## Other parameters
int.on = T #Intervention turned on at 2015
truncated.norm = F #whether to use truncated normal (instead of inverse gamma) as a proposal dist for gamma_b (hyperparameter for lambda_b)

## The values of the parameters the modelers gave us. True ??:
## see definitions in the draft
beta.m <- 0.6 # rate of transmission, between 0 - 1
InitialPropInfect.m <- 0.005 # the initial prevalence, bet 0 - 1
circatbirth.m <- 0.25 # c in the model, bet 0 - 1 used to specify initials and also proportion entering circumcised group at birth
epsilonart.m <- 0.7 # \epsilon_A, bet 0 - 1
epsiloncirc.m <- 0.6 # \epsilon_c, bet 0 - 1
mu.m <- 1.0 / 35 # \mu_1, can consider between 0 - 1
growth.m <- 0.03 # \epsilon, pop growth can consider bet 0 - 1
prograte.m <- 0.2 # \sigma_1, bet 0 - 1
infdeath.m <- 0.2  # \mu_2, bet 0 - 1
ARTdeath.m <- 0.05 # \mu_3, bet 0 - 1

ini.parms <- c(beta.m, InitialPropInfect.m, circatbirth.m) 
fixed.parmsfit <- c(epsilonart.m, epsiloncirc.m, mu.m, growth.m, prograte.m, infdeath.m, ARTdeath.m)
u.m <- c(fixed.parmsfit, ini.parms) ## the full set of parameters 

## note that \tau is optimized using other parameters through the fitmodel function for each location


#### 2. GENERATE FIELD DATA 

Z.input <- cbind(prev13, matrix(rep(u.m, each = length(prev13)), nrow = length(prev13))) ## appending the field data X with calibration parameters u
dat.true <- data.generate(nloc, prev13, u.start, time.list, 
                          Interventioninput, sample.size, bias.pct) ## the field data, the model output and the bias at field data augmented with true u


#### 3. GENERATING ONE INSTANCE OF MODEL OUTPUTS TO FIT GAUSSIAN PROCESS PARAMETERS

## getting the prior information
prior.u = get(prior)

### generate a parameter set from prior
u.n = prior.u(u.m, delta.prec, delta)

## Model output at communities for the newly generated u
mod.data.int <- mod.output(nloc, prev13, u.n, time.list, Interventioninput, int.on = T)
yM_int.n = drop(mod.data.int$yM_star)
yM_bas.n = drop(mod.data.bas$yM_star)
b_int.n = yF_int - yM_int.n
b_bas.n = yF_bas - yM_bas.n

X.new<-as.matrix(prev13)
GP_int.x = GP_fit(X.new, b_int.n, nug_thres = 10, corr = list(type = 'exponential', power = 2))
GP_bas.x = GP_fit(X.new, b_bas.n, nug_thres = 10, corr = list(type = 'exponential', power = 2))
beta_int.x = GP_int.x$beta
beta_bas.x = GP_bas.x$beta
mu_int.x = predict.GP(GP_int.x, X.new)$Y_hat
mu_bas.x = predict.GP(GP_bas.x, X.new)$Y_hat
R_int.x = corr_matrix(X.new, GP_int.x$beta, corr = list(type = 'exponential', power = 2))
R_bas.x = corr_matrix(X.new, GP_bas.x$beta, corr = list(type = 'exponential', power = 2))
gamma_int.x = 5 * 1 / GP_int.x$sig2 ## setting gamma for prior specification to = 5*lambdab_hat
gamma_bas.x = 5 * 1 / GP_bas.x$sig2 ## setting gamma for prior specification to = 5*lambdab_hat
Sig_int.x = GP_int.x$sig2 * R_int.x
Sig_bas.x = GP_bas.x$sig2 * R_bas.x

if (int.on == T){
  gammab = gamma_int.x
  betab = beta_int.x
  mub = rep(0, nloc)
  yF = yF_int
  R = corr_matrix(prev13, betab, corr = list(type = 'exponential', power = 2))
  R.fixed = R
} else {
  gammab = gamma_bas.x
  betab = beta_bas.x
  mub = rep(0, nloc)
  yF = yF_bas
  R.fixed = NULL
}


#### 4. RUNNING THE MCMC

c_beta2 = 1.5 * c_beta
c_beta_n = rep(c_beta, length(u.m))
c_beta_n[c(3,4,7)] = c_beta2

mcmc.run <- MCMC.ana.mh(prev13 = prev13, yF = yF, mod.output = mod.output, u.m = u.m, mub = mub, gammab = gammab,
                      betab = betab, time.list = time.list, sample.size = sample.size, scale1 = scale1,
                      scale2 = scale2, c_beta = c_beta_n, c_invgam = c_invgam, truncated.norm = truncated.norm,
                      sequential = T, delta = delta, delta.prec = delta.prec, int.on = int.on,
                      Interventioninput = Interventioninput, MCMC.samples = MCMC.samples, burn.in = burn.in,
                      thin.in = thin.in, update.theta = update.theta, prior = prior)









