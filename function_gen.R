####################################################################################
#### CODES FOR DIFFERENT FUNCTIONS TO RUN THE BAYESIAN VALIDATION FRAMEWORK
####################################################################################

## 1. TRUNCATED LOGIT TRANSFORMATION WITH A GIVEN PRECISION

logit.trunc <- function(p, prec.logit) {
  # boundary conditions
  p <- ifelse(p == 0, 0.5 * prec.logit, ifelse(p == 1, 1 - 0.5 * prec.logit, p))
  logit <- log(p) - log(1 - p)
  return(logit)
}

## 2. INVERSE LOGISTIC FUNCTION

expit <- function(x) 1 / (1 + exp( - x))


## 3. GENERATE FROM LOGITGAUSS PRIOR

logitgauss.prior <- function(u, prec = 1e - 2, delta = NULL){
  
  # slightly different way we set up variance for u from the uniform formulation
  if (is.null(delta) == T) delta <- prec * u * (1 - u)
  # we use the logit transform as most parameters are rates between 0 and 1
  logit.u <- log(u) - log(1 - u)
  # get transformed variance
  var.logit <- delta / ((u ^ 2) * ((1 - u) ^ 2))
  draw.logitu <- rnorm(length(u), logit.u, var.logit)
  urand <- expit(draw.logitu)
  urand <- ifelse(urand == 0, 1e - 5, ifelse(urand == 1, 1 - 1e - 5, urand))
  return(urand)
}

## 4. GENERATE FROM UNIFORM PRIOR

unif.prior <- function(u, prec = 1e - 2, delta = NULL){
  # the window depends on how close to 0 or 1 the mean is
  if (is.null(delta) == T) delta <- prec * (u * (1 - u) * (u < 0.9 & u > 0.1) + pmin(1 - u, u) * (u > 0.9 | u < 0.1))
  urand <- runif(length(u), u - delta, u + delta) 
  return(urand)
}

## 5. CALCULATE LAMBDA_F INVERSE

p.lambdaf_inv <- function(ym, b, sample.size = 250) 1 / (sample.size * expit(ym + b) * (1 - expit(ym + b)))

## 6. GENERATE LAMBDA_B INVERSE FROM GAMMA PRIOR

p.lambdab_inv <- function(gammab, n, b, R, prior = T){
  if (prior == T) lambdab <- rgamma(1, 1, scale = gammab)
  if (prior == F) lambdab<-rgamma(1, 1 + n / 2, rate = 1 / gammab + b %*% solve(R, t(b)))
  return(1 / lambdab)
}

## 7. FUNCTION TO GENERATE FIELD DATA (WITH/WITHOUT BIAS IN MODEL OUTPUTS)

data.generate <- function(nloc, prev13, u, time.list, Interventioninput, sample.size, bias.pct){
    
  bias.pct = as.numeric(strsplit(bias.pct, "[.]")[[1]])
  
  ## Define model parameters
  fixed.parmsfit <- u[1:7]
  beta <- u[8]
  InitialPropInfect <- u[9]
  circatbirth <- u[10]
  
  ## Running the model to get optimum parameters for 2013
  parameters <- fitmodel.n(nloc, prev13, fixed.parmsfit, beta, circatbirth, InitialPropInfect) #fit the model 
  
  InitialPopSizeadjusted <- rep(1, nloc)
  Popsizestartyearint <- rep(1, nloc)
  Adjustmentfactor <- rep(1, nloc)
  times <- seq(0, MaxTime, length = ((MaxTime * Nostepsinyear) + 1))  #time series for the output- 50 years, time step 0.05
  
  ## Rerun and get the population size in the start year of the intervention
  for (c in 1:nloc){
    parmsfit <- c(notatrisk = parameters[1, c], ART = 0, scale = 1, circmax = circatbirth, inttime = (MaxTime + 10), propcircbirth = circatbirth,
                  betafit = beta, trend = trenduncertain) 
    xstart <- c(S1 = (parameters[1, c] * circatbirth * (InitialPopSize * (1 - InitialPropInfect))), 
                S2 = ((1 - parameters[1, c]) * circatbirth * (InitialPopSize * (1 - InitialPropInfect))), 
                S3 = ((1 - parameters[1, c]) * (1 - circatbirth) * (InitialPopSize * (1 - InitialPropInfect))),
                S4 = (parameters[1, c] * (1 - circatbirth) * (InitialPopSize * (1 - InitialPropInfect))),
                I1 = (InitialPopSize * InitialPropInfect * 0.5), I2 = (InitialPopSize * InitialPropInfect * 0.5), AI = 0, NAI = 0) 
    output <- sim.mod.n(xstart, parmsfit, fixed.parmsfit, times) 
    
    #sum the number of individuals in each infection state in the year 2015
    Popsizestartyearint[c] <- sum(output[701, 2:9]) #year 35=output point 701 start year of int
    
    #making them all same as the first location.
    Adjustmentfactor[c] <- Popsizestartyearint[1] / Popsizestartyearint[c]
    InitialPopSizeadjusted[c] <- InitialPopSize * Adjustmentfactor[c]
  }
  
  ## Final results from the model
  final <- getresults.simu(parameters, Interventioninput, InitialPopSizeadjusted, fixed.parmsfit, beta, circatbirth, InitialPropInfect) 
  
  #deriving prevalence at different time points
  y_hat <- matrix(0, nloc, length(time.list))
  y_hat_b <- matrix(0, nloc, length(time.list))
  y_bins <- array(0, c(8, nloc, length(time.list)))
  y_bins_b <- array(0, c(8, nloc, length(time.list)))
  y_pred <- matrix(0, nloc, length(time.list))
  y_pred_b <- matrix(0, nloc, length(time.list))
  for (j in 1:length(time.list)){
    y_hat[,j] <- colSums(final$interventions[time.list[j], 6:7, 1:6]) / colSums(final$interventions[time.list[j], 2:9, 1:6])
    y_hat_b[,j] <- colSums(final$baseline[time.list[j], 6:7, 1:6]) / colSums(final$baseline[time.list[j], 2:9, 1:6])
  }
  
  # add bias at field observations 
  yM_int <- y_hat
  yM_bas <- y_hat_b
  
  yM_star_int <- logit.trunc(yM_int, prec.logit = 1 / sample.size)
  yM_star_bas <- logit.trunc(yM_bas, prec.logit = 1 / sample.size)
  
  yR_int <- pmax(pmin(y_hat + y_hat * (bias.pct / 100), 1), 0)
  yR_bas <- pmax(pmin(y_hat_b + y_hat_b * (bias.pct / 100), 1) ,0)
  yR_star_int <- logit.trunc(yR_int, prec.logit = 1 / sample.size)
  yR_star_bas <- logit.trunc(yR_bas, prec.logit = 1 / sample.size)
  
  b_int <- yR_star_int - yM_star_int
  b_bas <- yR_star_bas - yM_star_bas
  
  for (j in 1:length(time.list)){
    y_pred[, j] <- rbinom(nloc, sample.size, yR_int[,j]) / sample.size
    y_pred_b[, j] <- rbinom(nloc,sample.size,yR_bas[,j]) / sample.size
  }
  
  ## Field data
  yF_int <- y_pred
  yF_bas <- y_pred_b
  
  yF_star_int <- logit.trunc(y_pred, prec.logit = 1 / sample.size)
  yF_star_bas <- logit.trunc(y_pred_b, prec.logit = 1 / sample.size)
  
  
  return(list(yR_int = yR_int, yR_bas = yR_bas, yM_int = yM_int, yM_bas = yM_bas, yF_int = yF_int, 
              yF_bas = yF_bas, yR_star_int = yR_star_int, yR_star_bas = yR_star_bas, yM_star_int = yM_star_int, 
              yM_star_bas = yM_star_bas, yF_star_int = yF_star_int, yF_star_bas = yF_star_bas, b_int = b_int, 
              b_bas = b_bas, pars.ini = parameters))
  
}


## 8. MODEL OUTPUT FOR A GIVEN SET OF PARAMETERS

mod.output <- function(nloc, prev13, u, time.list, Interventioninput, int.on){
  
  fixed.parmsfit <- u[1:7]
  beta <- u[8]
  InitialPropInfect <- u[9]
  circatbirth <- u[10]
  
  ## Running the model to get optimum parameters for 2013
  parameters <- fitmodel.n(nloc, prev13, fixed.parmsfit,beta, circatbirth, InitialPropInfect) #fit the model 
  
  InitialPopSizeadjusted <- rep(1, nloc)
  Popsizestartyearint <- rep(1, nloc)
  Adjustmentfactor <- rep(1, nloc)
  times <- seq(0, MaxTime, length = ((MaxTime * Nostepsinyear) + 1))  #time series for the output- 50 years, time step 0.05
  
  #rerun and get the population size in the start year of the intervention
  for(c in 1:nloc){
    parmsfit <- c(notatrisk = parameters[1, c], ART = 0, scale = 1, circmax = circatbirth, inttime = (MaxTime + 10), propcircbirth = circatbirth,
                  betafit = beta, trend = trenduncertain)
    xstart <- c(S1 = (parameters[1, c] * circatbirth * (InitialPopSize * (1 - InitialPropInfect))), 
                S2 = ((1 - parameters[1, c]) * circatbirth * (InitialPopSize * (1 - InitialPropInfect))),
                S3 = ((1 - parameters[1, c])*(1 - circatbirth) * (InitialPopSize * (1 - InitialPropInfect))),
                S4 = (parameters[1, c] * (1 - circatbirth) * (InitialPopSize * (1 - InitialPropInfect))),
                I1 = (InitialPopSize * InitialPropInfect * 0.5), I2 = (InitialPopSize * InitialPropInfect * 0.5), AI = 0,NAI = 0) 
    output <- sim.mod.n(xstart, parmsfit, fixed.parmsfit, times) 
    #sum the number of individuals in each infection state in the year 2015
    Popsizestartyearint[c] <- sum(output[701, 2:9]) #year 35=output point 701 start year of int
    #making them all same as the first location.
    Adjustmentfactor[c] <- Popsizestartyearint[1] / Popsizestartyearint[c]
    InitialPopSizeadjusted[c] <- InitialPopSize * Adjustmentfactor[c]
  }
  
  ## Final results from the model
  final <- getresults.simu(parameters, Interventioninput, InitialPopSizeadjusted, fixed.parmsfit, beta, circatbirth, InitialPropInfect) #Output Cost/QALY time series
  
  if (int.on == T){
    y_hat <- matrix(0, nloc, length(time.list))
    for (j in 1:length(time.list)){
      y_hat[,j] <- colSums(final$interventions[time.list[j], 6:7, 1:6]) / colSums(final$interventions[time.list[j], 2:9, 1:6])
    }
    
  } else {
    y_hat <- matrix(0, nloc, length(time.list))
    for (j in 1:length(time.list)){
      y_hat[, j] <- colSums(final$baseline[time.list[j], 6:7, 1:6]) / colSums(final$baseline[time.list[j], 2:9, 1:6])
    }
    
  }
  y_hat_star <- logit.trunc(y_hat, prec.logit = 1 / sample.size)
  
  return(list(yM = y_hat, yM_star = y_hat_star))
  
}


## 9. MCMC ANALYSIS


MCMC.ana.mh <- function(prev13, yF, mod.output, u.m, mub, gammab, betab, time.list, sample.size, scale1, scale2, c_beta, c_invgam = NULL, 
                      truncated.norm = T, sequential = T, delta = NULL, delta.prec = 1e-2, int.on = T, Interventioninput = c(0.10, 0.05, 0.15), 
                      MCMC.samples = 5000, burn.in = 2000, thin.in = 1, update.theta = T, prior = 'logitgauss.prior'){
  
  ###############################
  ### step 0: Initial setup  ####
  ###############################
  
  ### 0. Get the prior
  prior.u <- get(prior)
  nloc <- length(prev13)
  
  ### 1. draw U from p(U)
  u.0 <- prior.u(u.m, delta.prec, delta)
  
  ### 2. model output
  mod.data.0 <- mod.output(nloc, prev13, u.0, time.list, Interventioninput, int.on)
  yM.0 <- drop(mod.data.0$yM_star)
  
  ### 3. lambda_b
  lb.inv.0 <- p.lambdab_inv(gammab, n = length(yF), b = NULL, R = NULL, prior = T)
  
  ### covariance matrix for bias b
  R <- corr_matrix(prev13, betab,  corr = list(type = 'exponential', power = 2))
  R.fixed = R

  Sig_b <- lb.inv.0*R

  #### 4. bias at U_0
  b.0 <- rmvnorm(1, mub, Sig_b)
  
  ### 3. Initial value for lambda_f
  lf.inv.0 <- p.lambdaf_inv(yM.0, b.0, sample.size)
  
  
  ###############################
  ### step i: MCMC analysis  ####
  ###############################  
  
  #inital values
  u.ini <- u.0
  yM.ini <- yM.0
  b.ini <- b.0
  lf.inv.ini <- lf.inv.0
  lb.inv.ini <- lb.inv.0
  
  par.ini <- c(u.ini, b.ini, lb.inv.ini)
  
  #the MCMC analysis
  if (update.theta == T) {
    mcmc.ana <- metro.hast(theta = par.ini, lik.fun = loglike.data, prior.fun = loglike.prior, scale1 = scale1, scale2 = scale2, c_beta = c_beta, u.m = u.m, delta.prec = delta.prec, delta = delta, c_invgam = c_invgam, 
                         truncated.norm = truncated.norm, mcmc.sample = MCMC.samples, burn.in = burn.in, thin.in = thin.in, sequential = sequential, prior = prior, R.fixed = R.fixed, gammab = gammab, mub = mub, betab = betab, 
                         sample.size = sample.size, nloc = nloc, prev13 = prev13, time.list = time.list, Interventioninput = Interventioninput, int.on = int.on)
  } else {
    mcmc.ana <- metro.hast.notheta(theta = par.ini, lik.fun = loglike.data, prior.fun = loglike.prior, scale1 = scale1, u.m = u.m, delta.prec = delta.prec, delta = delta, c_invgam = c_invgam, 
                                 mcmc.sample = MCMC.samples, burn.in = burn.in, thin.in = thin.in, sequential = sequential, prior = prior, R.fixed = R.fixed, gammab = gammab, mub = mub, betab = betab, 
                                 sample.size = sample.size, nloc = nloc, prev13 = prev13, time.list = time.list, Interventioninput = Interventioninput, int.on = int.on)
  }
  
  return(mcmc.ana)
}


## 10. CALCULATE LOG LIKELIHOOD FOR THE DATA GIVEN MODEL PARAMETERS

loglike.data <- function(par){
  
  library(mvtnorm)
  u <- par[1:length(u.m)]
  b <- par[(length(u.m) + 1):(length(u.m) + length(prev13))]
  
  ### model output at u
  mod.data <- mod.output(nloc, prev13, u, time.list, Interventioninput, int.on)
  yM <- drop(mod.data$yM_star)
  
  ### covariance matrix for yF
  lf.inv <- p.lambdaf_inv(yM, b, sample.size)
  Sig_yF <- diag(as.vector(lf.inv))
  
  ### likelihood
  l <- log(dmvnorm(yF, t(yM) + b, Sig_yF))
  return(list(loglike = l,yM = yM))
}


## 11. CALCULATE LOG LIKELIHOOD FOR THE MODEL PARAMETERS


loglike.prior <- function(par){
  
  library(mvtnorm)
  u <- par[1:length(u.m)]
  b <- par[(length(u.m) + 1):(length(u.m) + length(prev13))]
  lb.inv <- par[length(par)]
  
  ### covariance matrix for bias b
  R.new <- R.fixed
  if (is.null(R.new) == T) {
    Z = cbind(prev13, matrix(rep(u, each = length(prev13)), nrow = length(prev13)))
    R.new <- corr_matrix(Z, betab,  corr = list(type = 'exponential', power = 2))
  }
  Sig_b <- lb.inv * R.new # covariance matrix for bias b
  
  ### likelihood
  if (prior == 'unif.prior'){
    if (is.null(delta) == T) delta <- u.m * (1 - u.m) * (u.m < 0.9 & u.m > 0.1)  +  pmin(1 - u.m, u.m)*(u.m > 0.9 | u.m < 0.1)
    l1 <- sum(log(dunif(u, u.m - delta.prec * delta, u.m + delta.prec * delta)))
  } else if(prior == 'logitgauss.prior'){
    if (is.null(delta) == T) delta <- delta.prec * u.m * (1 - u.m)
    logit.u <- log(u) - log(1 - u)
    mean.logit.u <- log(u.m) - log(1 - u.m)
    var.logit.u <- delta / ((u.m ^ 2) * ((1 - u.m)^2))
    l1 <- sum(log(dnorm(logit.u, mean.logit.u, var.logit.u)))
  } 
  l2 <- log(dmvnorm(b, mub, Sig_b))
  l3 <- log(dgamma(1 / lb.inv, 1, scale = gammab))
  
  return(l1 + l2 + l3)
}


## 12. METROPOLIS HASTINGS ALGORITHM WITH ALL PARAMETERS BEING UPDATED SIMULTANEOUSLY

metro.hast <- function(theta, lik.fun, prior.fun, scale1, scale2, c_beta, u.m, delta.prec, delta, c_invgam, truncated.norm, mcmc.sample, burn.in, 
                       thin.in, sequential, prior, R.fixed, gammab, mub, betab, sample.size, nloc, prev13, time.list, Interventioninput, int.on) {
  prior.u <- get(prior)
  if (! sequential) {
    if (! require(mvtnorm)) stop("You need to install the mvtnorm library.")
  }
  if (! truncated.norm & is.null(c_invgam)) stop("need to specify c_invgam for inverse gamma proposal")
  if (length(c_beta) == 1) c_beta = rep(c_beta, length(u.m))
  if (length(scale1) == 1) scale1 = rep(scale1, length(prev13))
  eff.samp <- burn.in + mcmc.sample
  ndim <- length(theta)
  chain <- matrix(NA, nrow = eff.samp, ncol = ndim)
  modout<- matrix(0, nrow = eff.samp, ncol = length(prev13))
  lik <- vector(length = eff.samp)
  last <- theta
  #print(last)
  chain[1, ] <- theta
  lik.summ <- lik.fun(theta)
  last.lik <- lik.summ$loglike
  lik[1] <- last.lik
  modout[1, ] <- lik.summ$yM
  last.prior <- prior.fun(theta)
  it <- 1
  if (is.finite(last.lik)) {
    for (it in seq(2, eff.samp)) {
      print(it)
      if (sequential) {
        for (i in 1:length(theta)) {
          accept <- FALSE
          proposal <- last
          proposal[i] <- proposal.uni(i, last[i], scale1, c_beta, scale2, c_invgam, truncated.norm)
          qratio.proposal <- proposal.loglike.uni(i, last[i], proposal[i], scale1, c_beta, scale2, c_invgam, truncated.norm)
          # print(qratio.proposal)
          proposal.prior <- prior.fun(proposal)
          # print(proposal.prior)
          if (is.finite(proposal.prior)) {
            proposal.summ <- lik.fun(proposal)
            proposal.lik <- proposal.summ$loglike
            alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
            if (alpha > runif(1)) accept <- TRUE
          }
          #print(proposal)
          #print(accept)
          if (accept) {
            last <- proposal
            lik.summ <- proposal.summ
            last.lik <- lik.summ$loglike
            last.prior <- proposal.prior
          }
          #print(last)
          #print(i)
        }
        chain[it, ] <- last
        lik[it] <- last.lik
        modout[it, ] <- lik.summ$yM
        
      } else {
        
        accept <- FALSE
        proposal <- proposal.mult(last, scale1, c_beta, scale2, c_invgam, truncated.norm)
        qratio.proposal <- proposal.loglike.mult(last, proposal, scale1, c_beta, scale2, c_invgam, truncated.norm)
        proposal.prior <- prior.fun(proposal)
        if (is.finite(proposal.prior)) {
          proposal.summ <- lik.fun(proposal)
          proposal.lik <- proposal.summ$loglike
          alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
          if (alpha > runif(1)) accept <- TRUE
        }
        if (accept) {
          last <- chain[it, ] <- proposal
          last.lik <- lik[it] <- proposal.lik
          modout[it, ] <- proposal.summ$yM
          last.prior <- proposal.prior
        } else {
          last <- chain[it, ] <- chain[it-1, ]
          last.lik <- lik[it] <- lik[it-1]
          modout[it, ] <- modout[it-1, ]
        }
      }
    }
  }
  mcmcid <- seq(1, MCMC.samples, thin.in) + burn.in
  accept.vector <- apply(chain, 2, function(x) length(unique(x)) / length(x))
  mcmc.chain <- chain[mcmcid, ]
  return(list(full.chain = chain, lik = lik, mcmc.chain = mcmc.chain, model.output = modout, accept.vector = accept.vector))
}


## 13. METROPOLIS HASTINGS ALGORITHM WITH NO PARAMETERS BEING UPDATED (NOTE THAT PRIOR FOR DISCREPANCY IS STILL UPDATED)

metro.hast.notheta <- function(theta, lik.fun, prior.fun, scale1, u.m, delta.prec, delta, c_invgam, mcmc.sample, burn.in, thin.in, sequential, 
                               prior, R.fixed, gammab, mub, betab, sample.size, nloc, prev13, time.list, Interventioninput, int.on) {
  prior.u<-get(prior)
  if (! sequential) {
    if (! require(mvtnorm)) stop("You need to install the mvtnorm library.")
  }
  if (length(scale1) == 1) scale1 = rep(scale1, length(prev13))
  eff.samp<-burn.in + mcmc.sample
  ndim <- length(theta)
  chain <- matrix(NA, nrow = eff.samp, ncol = ndim)
  modout<- matrix(0, nrow = eff.samp, ncol = length(prev13))
  lik <- vector(length = eff.samp)
  last <- theta
  #print(last)
  chain[1, ] <- last
  last.summ <- lik.fun(last)
  last.lik <- last.summ$loglike
  modout[1, ]<-last.summ$yM
  lik[1] <- last.lik
  last.prior <- prior.fun(last)
  it <- 1
  if (is.finite(last.lik)) {
    for (it in seq(2, eff.samp)) {
      print(it)
      if (sequential) {
        last[1:length(u.m)] <- prior.u(u.m, delta.prec, delta)
        last.prior <- prior.fun(last)
        last.summ <- lik.fun(last)
        last.lik <- last.summ$loglike
        for (i in (length(u.m) + 1):length(theta)) {
          accept <- FALSE
          proposal <- last
          proposal[i] <- proposal.uni.notheta(i, last[i], scale1, c_invgam, u.m, prev13)
          qratio.proposal <- proposal.loglike.uni.notheta(i, last[i], proposal[i], scale1, c_invgam, u.m, prev13)
          # print(qratio.proposal)
          proposal.prior <- prior.fun(proposal)
          # print(proposal.prior)
          if (is.finite(last.prior) & is.finite(proposal.prior)) {
            proposal.summ <- lik.fun(proposal)
            proposal.lik <- proposal.summ$loglike
            if (is.infinite(last.lik) & is.infinite(proposal.lik) & last.lik == proposal.lik){
              alpha <- exp((proposal.prior - last.prior) + qratio.proposal)
            } else{
              alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
            }
            if (alpha > runif(1)) accept <- TRUE
          }
          #print(proposal)
          #print(accept)
          if (accept) {
            last <- proposal
            last.summ <- proposal.summ
            last.lik <- last.summ$loglike
            last.prior <- proposal.prior
          }
          #print(last)
          #print(i)
        }
        chain[it, ] <- last
        lik[it] <- last.lik
        modout[it, ] <- last.summ$yM
        
      } else {
        last[1:length(u.m)] <- prior.u(u.m, delta.prec, delta)
        last.prior <- prior.fun(last)
        last.summ <- lik.fun(last)
        last.lik <- last.summ$loglike
        last.u <- last[1:length(u.m)]
        last.nou <- last[(length(u.m) + 1):length(last)]
        accept <- FALSE
        proposal.nou <- proposal.mult.notheta(last.nou, scale1, c_invgam, u.m, prev13)
        proposal <- c(last.u, proposal.nou)
        qratio.proposal <- proposal.loglike.mult.notheta(last.nou, proposal.nou, scale1, c_invgam, u.m, prev13)
        proposal.prior <- prior.fun(proposal)
        if (is.finite(proposal.prior)) {
          proposal.summ <- lik.fun(proposal)
          proposal.lik<- proposal.summ$loglike
          if (is.infinite(last.prior) & is.infinite(proposal.prior) & last.prior == proposal.prior){
            alpha <- exp((proposal.prior - last.prior) + qratio.proposal)
          } else{
            alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
          }
          if (alpha > runif(1)) accept <- TRUE
        }
        if (accept) {
          last <- chain[it, ] <- proposal
          last.lik <- lik[it] <- proposal.lik
          modout[it, ] <- proposal.summ$yM
          last.prior <- proposal.prior
        } else {
          chain[it, ] <- last
          lik[it] <- last.lik
          modout[it, ] <- last.summ$yM
          
        }
      }
    }
  }
  mcmcid <- seq(1, MCMC.samples, thin.in) + burn.in
  accept.vector <- apply(chain, 2, function(x) length(unique(x)) / length(x))
  mcmc.chain <- chain[mcmcid, ]
  return(list(full.chain = chain, lik = lik, mcmc.chain = mcmc.chain, model.output = modout, accept.vector = accept.vector))
}

## 14. GENERATE FROM PROPOSAL DISTRIBUTION WHEN ALL PARAMETERS ARE BEING UPDATED SIMULTANEOUSLY (STEP BY STEP UNIVARIATE)

proposal.uni <- function(i, theta1, scale1, c_beta, scale2, c_invgam, truncated.norm){
  if (i %in% 1:length(u.m)){
    repeat{
      proposal <- rbeta(1, shape1 = c_beta[i] * theta1, shape2 = c_beta[i] * (1 - theta1))
      if(proposal < 1 | proposal > 0)
        break
    }
  } else if (i %in% (length(u.m) + 1):(length(u.m) + length(prev13))){
    proposal <- rnorm(1, theta1, sqrt(scale1[i - length(u.m)]))
    
  } else if (i %in% (length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(gammab))){
    if (truncated.norm) {
      repeat{
        proposal <- rnorm(1, theta1, sqrt(scale2))
        if (proposal > 0) 
          break
      }
    } else {
      proposal <- rinvgamma(1, shape = c_invgam + 1, scale = c_invgam * theta1)
    }
  }
  return(proposal)
}

## 15. GENERATE FROM PROPOSAL DISTRIBUTION WHEN NO PARAMETERS ARE BEING UPDATED (NOTE THAT PRIOR FOR DISCREPANCY IS STILL UPDATED) (STEP BY STEP UNIVARIATE)

proposal.uni.notheta <- function(i, theta1, scale1, c_invgam, u.m, prev13){
  if (i %in% (length(u.m) + 1):(length(u.m) + length(prev13))){
    proposal <- rnorm(1, theta1, sqrt(scale1[i - length(u.m)]))
    
  } else if (i %in% (length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(gammab))){
    proposal <- rinvgamma(1, shape = c_invgam + 1, scale = c_invgam * theta1)
  }
  return(proposal)
}

## 16. CALCULATE Q RATIO (LOG LIKELIHOOD RATIO) FOR THE PROPOSAL VALUE TO THE PREVIOUS VALUE (STEP BY STEP UNIVARIATE) WHEN ALL PARAMETERS ARE BEING UPDATED

proposal.loglike.uni <- function(i, theta1, theta2, scale1, c_beta, scale2, c_invgam, truncated.norm){
  if (i %in% 1:length(u.m)){
    qratio.proposal  <-  dbeta(theta1, shape1 = c_beta[i] * theta2, shape2 = c_beta[i] * (1 - theta2), log = T) -
      dbeta(theta2, shape1 = c_beta[i] * theta1, shape2 = c_beta[i]*(1-theta1), log = T)
  } else if (i %in% (length(u.m) + 1):(length(u.m) + length(prev13))){
    qratio.proposal  <-  0
  } else if (i %in% (length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(gammab))) {
    if (truncated.norm) {
      qratio.proposal <- pnorm(theta1, mean = 0, sd = sqrt(scale2), log.p = T) - pnorm(theta2, mean = 0, sd = sqrt(scale2), log.p = T)
    } else {
      qratio.proposal <- log(dinvgamma(theta1, shape = c_invgam + 1, scale = c_invgam * theta2)) - 
        log(dinvgamma(theta2, shape = c_invgam + 1, scale = c_invgam * theta1))
    }
    
    
  }
  return(qratio.proposal)
  
}

## 17. CALCULATE Q RATIO (LOG LIKELIHOOD RATIO) FOR THE PROPOSAL VALUE TO THE PREVIOUS VALUE (STEP BY STEP UNIVARIATE) WHEN NO PARAMETERS ARE BEING UPDATED

proposal.loglike.uni.notheta <- function(i, theta1, theta2, scale1, c_invgam, u.m, prev13){
  if (i %in% (length(u.m) + 1):(length(u.m) + length(prev13))){
    qratio.proposal  <-  0
  } else if (i %in% (length(u.m) + length(prev13)):(length(u.m) + length(prev13) + length(gammab))) {
    qratio.proposal <- log(dinvgamma(theta1, shape = c_invgam + 1, scale = c_invgam * theta2)) - 
      log(dinvgamma(theta2, shape = c_invgam + 1, scale = c_invgam * theta1))
  }
  return(qratio.proposal)
  
}

## 18. GENERATE FROM PROPOSAL DISTRIBUTION WHEN ALL PARAMETERS ARE BEING UPDATED SIMULTANEOUSLY (TOGETHER MULTIVARIATE)

proposal.mult <- function(theta, scale1, c_beta, scale2, c_invgam, truncated.norm){
  
  theta.u <- theta[1:length(u.m)]
  theta.u.n <- rep(0, length(theta.u))
  for (i in 1:length(theta.u)){
    repeat{
      theta.u.n[i] = rbeta(1, shape1 = c_beta[i] * theta.u[i], shape2 = c_beta[i] * (1-theta.u[i]))
      if(theta.u.n[i] < 1 | theta.u.n[i] > 0)
        break
    }
  }
  theta.b <- theta[(length(u.m) + 1):(length(u.m) + length(prev13))]
  V1 <- diag(scale1)
  theta.b.n <- rmvnorm(1, mean = theta.b, sigma = V1)
  
  theta.lb <- theta[(length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(c_invgam))]
  if (truncated.norm) {
    repeat{
      theta.lb.n <- rnorm(1, theta.lb, sqrt(scale2))
      if (theta.lb.n > 0) 
        break
    }
  } else {
    theta.lb.n <- rinvgamma(1, shape = c_invgam + 1, scale = c_invgam * theta.lb)
  }
  
  return(c(theta.u.n, theta.b.n, theta.lb.n))  
}

## 19. CALCULATE Q RATIO (LOG LIKELIHOOD RATIO) FOR THE PROPOSAL VALUE TO THE PREVIOUS VALUE (TOGETHER MULTIVARIATE) WHEN ALL PARAMETERS ARE BEING UPDATED


proposal.loglike.mult <- function(theta1, theta2, scale1, c_beta, scale2, c_invgam, truncated.norm){
  
  theta1.u <- theta1[1:length(u.m)]
  theta2.u <- theta2[1:length(u.m)]
  loglike.u <- dbeta(theta1.u, shape1 = c_beta * theta2.u, shape2 = c_beta * (1 - theta2.u), log = T) -
    dbeta(theta2.u, shape1 = c_beta * theta1.u, shape2 = c_beta * (1 - theta1.u), log = T)
  
  loglike.b <- 0
  
  theta1.lb <- theta1[(length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(c_invgam))]
  theta2.lb <- theta2[(length(u.m) + length(prev13) + 1):(length(u.m) + length(prev13) + length(c_invgam))]
  loglike.lb <- log(dinvgamma(theta1.lb, shape = c_invgam + 1, scale = c_invgam*theta2.lb)) - 
    log(dinvgamma(theta2.lb, shape = c_invgam + 1, scale = c_invgam * theta1.lb))
  return(sum(loglike.u) + sum(loglike.b) + sum(loglike.lb))  
}


## 20. GENERATE FROM PROPOSAL DISTRIBUTION WHEN NO PARAMETERS ARE BEING UPDATED (NOTE THAT PRIOR FOR DISCREPANCY IS STILL UPDATED) (TOGETHER MULTIVARIATE)


proposal.mult.notheta <- function(theta, scale1, c_invgam, u.m, prev13){
  
  theta.b <- theta[1:length(prev13)]
  V1 <- diag(scale1)
  theta.b.n <- rmvnorm(1, mean = theta.b, sigma = V1)
  
  theta.lb <- theta[(length(prev13) + 1):(length(prev13) + length(gammab))]
  ### no need for truncated normal as we are not using that
  theta.lb.n <- rinvgamma(1, shape = c_invgam + 1, scale = c_invgam * theta.lb)
  
  return(c(theta.b.n, theta.lb.n))  
}

## 21. CALCULATE Q RATIO (LOG LIKELIHOOD RATIO) FOR THE PROPOSAL VALUE TO THE PREVIOUS VALUE (TOGETHER MULTIVARIATE) WHEN NO PARAMETERS ARE BEING UPDATED


proposal.loglike.mult.notheta <- function(theta1, theta2, scale1, c_invgam, u.m, prev13){
  
  loglike.b <- 0
  
  theta1.lb <- theta1[(length(prev13) + 1):(length(prev13) + length(gammab))]
  theta2.lb <- theta2[(length(prev13) + 1):(length(prev13) + length(gammab))]
  loglike.lb <- log(dinvgamma(theta1.lb, shape = c_invgam + 1, scale = c_invgam * theta2.lb)) - 
    log(dinvgamma(theta2.lb, shape = c_invgam + 1, scale = c_invgam * theta1.lb))
  return(sum(loglike.b) + sum(loglike.lb))  
}




####################################################################################################
## THE FUNCTIONS BELOW ARE FOR SETTING 2; THE ABOVE FUNCTIONS WERE USED IN SETTINGS 1 AND 3
####################################################################################################


## 22. MCMC ANALYSIS (SUBSET OF THETA)

MCMC.ana.mh.subtheta <- function(prev13, yF, mod.output, u.m, mub, gammab, betab, time.list, sample.size, scale1, scale2, 
                               c_beta, c_invgam = NULL, truncated.norm = T, sequential = T, delta = NULL, delta.prec = 1e-2, 
                               int.on = T, Interventioninput = c(0.10, 0.05, 0.15), MCMC.samples = 5000, burn.in = 2000, 
                               thin.in = 1, update.theta = T, prior = 'logitgauss.prior', thetaind = thetaind){
  
  ###############################
  ### step 0: Initial setup  ####
  ###############################
  
  ### 0. Get the prior
  prior.u <- get(prior)
  nloc <- length(prev13)
  
  likeval <- Inf
  while (!is.finite(likeval)){
    ### 1. draw U from p(U)
    u.m.ind <- u.m[thetaind]
    u.0 <- u.m
    u.0[thetaind] <- prior.u(u.m.ind, delta.prec, delta)
    
    ### 2. model output
    mod.data.0 <- mod.output(nloc, prev13, u.0, time.list, Interventioninput, int.on)
    yM.0 <- drop(mod.data.0$yM_star)
    
    ### 3. lambda_b
    lb.inv.0 <- p.lambdab_inv(gammab, n = length(yF), b = NULL, R = NULL, prior = T)
    
    ### covariance matrix for bias b
    R <- corr_matrix(prev13, betab,  corr = list(type = 'exponential', power = 2))
    R.fixed = R
    Sig_b <- lb.inv.0 * R
    
    
    #### 4. bias at U_0
    b.0 <- rmvnorm(1, mub, Sig_b)
    
    ### 3. Initial value for lambda_f
    lf.inv.0 <- p.lambdaf_inv(yM.0, b.0, sample.size)
    
    #inital values
    u.ini <- u.0
    yM.ini <- yM.0
    b.ini <- b.0
    lf.inv.ini <- lf.inv.0
    lb.inv.ini <- lb.inv.0
    
    par.ini <- c(u.ini[thetaind], b.ini, lb.inv.ini)
    
    likeval <- loglike.data.subtheta(par.ini)$loglike
    
  }
  
  #the MCMC analysis
  if (update.theta == T) {
    mcmc.ana <- metro.hast.subtheta(theta = par.ini, lik.fun = loglike.data.subtheta, prior.fun = loglike.prior.subtheta, scale1 = scale1, scale2 = scale2, c_beta = c_beta, u.m = u.m, delta.prec = delta.prec, delta = delta, c_invgam = c_invgam, 
                                  truncated.norm = truncated.norm, mcmc.sample = MCMC.samples, burn.in = burn.in, thin.in = thin.in, sequential = sequential, prior = prior, R.fixed = R.fixed, gammab = gammab, mub = mub, betab = betab, 
                                  sample.size = sample.size, nloc = nloc, prev13 = prev13, time.list = time.list, Interventioninput = Interventioninput, int.on = int.on, thetaind = thetaind)
  } else {
    mcmc.ana <- metro.hast.subtheta.noupdate(theta = par.ini, lik.fun = loglike.data.subtheta, prior.fun = loglike.prior.subtheta, scale1 = scale1, u.m = u.m, delta.prec = delta.prec, delta = delta, c_invgam = c_invgam, 
                                           mcmc.sample = MCMC.samples, burn.in = burn.in, thin.in = thin.in, sequential = sequential, prior = prior, R.fixed = R.fixed, gammab = gammab, mub = mub, betab = betab, 
                                           sample.size = sample.size, nloc = nloc, prev13 = prev13, time.list = time.list, Interventioninput = Interventioninput, int.on = int.on, thetaind = thetaind)
  }
  
  return(mcmc.ana)
}


## 23. METROPOLIS HASTINGS ALGORITHM WITH ALL PARAMETERS BEING UPDATED SIMULTANEOUSLY (SUBSET OF PARAMETERS)

metro.hast.subtheta <- function(theta, lik.fun, prior.fun, scale1, scale2, c_beta, u.m, delta.prec, delta, c_invgam, truncated.norm, mcmc.sample, burn.in, 
                                thin.in, sequential, prior, R.fixed, gammab, mub, betab, sample.size, nloc, prev13, time.list, Interventioninput, int.on, thetaind) {
  prior.u <- get(prior)
  if (! sequential) {
    if (! require(mvtnorm)) stop("You need to install the mvtnorm library.")
  }
  if (! truncated.norm & is.null(c_invgam)) stop("need to specify c_invgam for inverse gamma proposal")
  if (length(c_beta) == 1) c_beta = rep(c_beta, length(u.m[thetaind]))
  if (length(scale1) == 1) scale1 = rep(scale1, length(prev13))
  eff.samp <- burn.in + mcmc.sample
  ndim <- length(theta)
  chain <- matrix(NA, nrow = eff.samp, ncol = ndim)
  modout<- matrix(0, nrow = eff.samp, ncol = length(prev13))
  lik <- vector(length = eff.samp)
  last <- theta
  #print(last)
  chain[1, ] <- theta
  last.summ <- lik.fun(theta)
  last.lik <- last.summ$loglike
  lik[1] <- last.lik
  modout[1, ] <- last.summ$yM
  last.prior <- prior.fun(theta)
  it <- 1
  
  if (is.finite(last.lik)) {
    for (it in seq(2, eff.samp)) {
      print(it)
      if (sequential) {
        for (i in 1:length(theta)) {
          accept <- FALSE
          proposal <- last
          proposal[i] <- proposal.uni(i, last[i], scale1, c_beta, scale2, c_invgam, truncated.norm)
          qratio.proposal <- proposal.loglike.uni(i, last[i], proposal[i], scale1, c_beta, scale2, c_invgam, truncated.norm)
          # print(qratio.proposal)
          proposal.prior <- prior.fun(proposal)
          # print(proposal.prior)
          if (is.finite(proposal.prior)) {
            proposal.summ <- lik.fun(proposal)
            proposal.lik <- proposal.summ$loglike
            alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
            if (alpha > runif(1)) accept <- TRUE
          }
          #print(proposal)
          #print(accept)
          if (accept) {
            last <- proposal
            last.summ <- proposal.summ
            last.lik <- last.summ$loglike
            last.prior <- proposal.prior
          }
          #print(last)
          #print(i)
        }
        chain[it, ] <- last
        lik[it] <- last.lik
        modout[it, ] <- last.summ$yM
        
      } else {
        
        accept <- FALSE
        proposal <- proposal.mult(last, scale1, c_beta, scale2, c_invgam, truncated.norm)
        qratio.proposal <- proposal.loglike.mult(last, proposal, scale1, c_beta, scale2, c_invgam, truncated.norm)
        proposal.prior <- prior.fun(proposal)
        if (is.finite(proposal.prior)) {
          proposal.summ <- lik.fun(proposal)
          proposal.lik <- proposal.summ$loglike
          alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
          if (alpha > runif(1)) accept <- TRUE
        }
        if (accept) {
          last <- chain[it, ] <- proposal
          last.lik <- lik[it] <- proposal.lik
          modout[it, ] <- proposal.summ$yM
          last.prior <- proposal.prior
        } else {
          last <- chain[it, ] <- chain[it-1, ]
          last.lik <- lik[it] <- lik[it-1]
          modout[it, ] <- modout[it-1, ]
        }
      }
    }
  }
  mcmcid <- seq(1, MCMC.samples, thin.in) + burn.in
  accept.vector <- apply(chain, 2, function(x) length(unique(x)) / length(x))
  mcmc.chain <- chain[mcmcid, ]
  return(list(full.chain = chain, lik = lik, mcmc.chain = mcmc.chain, model.output = modout, accept.vector = accept.vector))
}

## 24. METROPOLIS HASTINGS ALGORITHM WITH ALL PARAMETERS BEING UPDATED SIMULTANEOUSLY (SUSBET OF PARAMETERS)

metro.hast.subtheta.noupdate <- function(theta, lik.fun, prior.fun, scale1, u.m, delta.prec, delta, c_invgam, mcmc.sample, burn.in, thin.in, sequential, prior, 
                                         R.fixed, gammab, mub, betab, sample.size, nloc, prev13, time.list, Interventioninput, int.on, thetaind) {
  prior.u <- get(prior)
  if (! sequential) {
    if (! require(mvtnorm)) stop("You need to install the mvtnorm library.")
  }
  if (length(scale1) == 1) scale1 = rep(scale1, length(prev13))
  eff.samp <- burn.in + mcmc.sample
  ndim <- length(theta)
  chain <- matrix(NA, nrow = eff.samp, ncol = ndim)
  modout <- matrix(0, nrow = eff.samp, ncol = length(prev13))
  lik <- vector(length = eff.samp)
  last <- theta
  #print(last)
  chain[1, ] <- last
  last.summ <- lik.fun(last)
  last.lik <- last.summ$loglike
  modout[1, ] <- last.summ$yM
  lik[1] <- last.lik
  last.prior <- prior.fun(last)
  it <- 1
  if (is.finite(last.lik)) {
    for (it in seq(2, eff.samp)) {
      print(it)
      if (sequential) {
        indlik = 1
        while (indlik == 1){
          last[1:length(thetaind)] <- prior.u(u.m[thetaind], delta.prec, delta)
          last.prior <- prior.fun(last)
          last.summ <- lik.fun(last)
          last.lik <- last.summ$loglike
          indlik = is.infinite(last.lik) 
        }
        for (i in (length(thetaind) + 1):length(theta)) {
          accept <- FALSE
          proposal <- last
          proposal[i] <- proposal.uni.notheta(i, last[i], scale1, c_invgam, u.m, prev13)
          qratio.proposal <- proposal.loglike.uni.notheta(i, last[i], proposal[i], scale1, c_invgam, u.m, prev13)
          # print(qratio.proposal)
          proposal.prior <- prior.fun(proposal)
          # print(proposal.prior)
          if (is.finite(last.prior) & is.finite(proposal.prior)) {
            proposal.summ <- lik.fun(proposal)
            proposal.lik <- proposal.summ$loglike
            alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
            if (alpha > runif(1)) accept <- TRUE
          }
          #print(proposal)
          #print(accept)
          if (accept) {
            last <- proposal
            last.summ <- proposal.summ
            last.lik <- last.summ$loglike
            last.prior <- proposal.prior
          }
          #print(last)
          #print(i)
        }
        chain[it, ] <- last
        lik[it] <- last.lik
        modout[it, ] <- last.summ$yM
        
      } else {
        indlik = 1
        while (indlik == 1){
          last[1:length(thetaind)] <- prior.u(u.m[thetaind], delta.prec, delta)
          last.prior <- prior.fun(last)
          last.summ <- lik.fun(last)
          last.lik <- last.summ$loglike
          indlik = is.infinite(last.lik) 
        }
        last.u <- last[1:length(thetaind)]
        last.nou <- last[(length(thetaind) + 1):length(last)]
        accept <- FALSE
        proposal.nou <- proposal.mult.notheta(last.nou, scale1, c_invgam, u.m, prev13)
        proposal <- c(last.u, proposal.nou)
        qratio.proposal <- proposal.loglike.mult.notheta(last.nou, proposal.nou, scale1, c_invgam, u.m, prev13)
        proposal.prior <- prior.fun(proposal)
        if (is.finite(proposal.prior)) {
          proposal.summ <- lik.fun(proposal)
          proposal.lik <- proposal.summ$loglike
          alpha <- exp((proposal.lik + proposal.prior - last.lik - last.prior) + qratio.proposal)
          if (alpha > runif(1)) accept <- TRUE
        }
        if (accept) {
          last <- chain[it, ] <- proposal
          last.lik <- lik[it] <- proposal.lik
          modout[it, ] <- proposal.summ$yM
          last.prior <- proposal.prior
        } else {
          chain[it, ] <- last
          lik[it] <- last.lik
          modout[it, ] <- last.summ$yM
          
        }
      }
    }
  }
  mcmcid <- seq(1, MCMC.samples, thin.in) + burn.in
  accept.vector <- apply(chain, 2, function(x) length(unique(x)) / length(x))
  mcmc.chain <- chain[mcmcid, ]
  return(list(full.chain = chain, lik = lik, mcmc.chain = mcmc.chain, model.output = modout, accept.vector = accept.vector))
}




## 25. CALCULATE LOG LIKELIHOOD FOR THE DATA GIVEN MODEL PARAMETERS (SUBSET OF PARAMETERS)

loglike.data.subtheta <- function(par){
  
  library(mvtnorm)
  u <- u.m
  u[thetaind] <- par[1:length(thetaind)]
  
  b <- par[(length(thetaind) + 1):(length(thetaind) + length(prev13))]
  
  ### model output at u
  mod.data <- mod.output(nloc, prev13, u, time.list, Interventioninput, int.on)
  yM <- drop(mod.data$yM_star)
  
  ### covariance matrix for yF
  lf.inv <- p.lambdaf_inv(yM, b, sample.size)
  Sig_yF <- diag(as.vector(lf.inv))
  
  ### likelihood
  l <- log(dmvnorm(yF, t(yM) + b , Sig_yF))
  return(list(loglike = l,yM = yM))
  
}

## 26. CALCULATE LOG LIKELIHOOD FOR THE MODEL PARAMETERS (SUBSET OF PARAMETERS)


loglike.prior.subtheta <- function(par){
  
  library(mvtnorm)
  u <- u.m
  u[thetaind] <- par[1:length(thetaind)]
  b <- par[(length(thetaind) + 1):(length(thetaind) + length(prev13))]
  lb.inv <- par[length(par)]
  
  ### covariance matrix for bias b
  R.new <- R.fixed
  if (is.null(R.new) == T) {
    Z = cbind(prev13, matrix(rep(u, each = length(prev13)), nrow = length(prev13)))
    R.new <- corr_matrix(Z, betab,  corr = list(type = 'exponential', power = 2))
  }
  Sig_b <- lb.inv * R.new # covariance matrix for bias b
  
  ### likelihood
  if (prior == 'unif.prior'){
    if (is.null(delta) == T) delta <- u.m[thetaind] * (1 - u.m[thetaind]) * (u.m[thetaind] < 0.9 & u.m[thetaind] > 0.1)  +  
        pmin(1 - u.m[thetaind], u.m[thetaind]) * (u.m[thetaind] > 0.9 | u.m[thetaind] < 0.1)
    l1 <- sum(log(dunif(u[thetaind], u.m[thetaind] - delta.prec * delta, u.m[thetaind] + delta.prec * delta)))
  } else if(prior == 'logitgauss.prior'){
    if (is.null(delta) == T) delta <- delta.prec * u.m[thetaind] * (1 - u.m[thetaind])
    logit.u <- log(u[thetaind]) - log(1 - u[thetaind])
    mean.logit.u <- log(u.m[thetaind]) - log(1 - u.m[thetaind])
    var.logit.u <- delta / ((u.m[thetaind] ^ 2) * ((1 - u.m[thetaind]) ^ 2))
    l1 <- sum(log(dnorm(logit.u, mean.logit.u, var.logit.u)))
  } 
  l2 <- log(dmvnorm(b, mub, Sig_b))
  l3 <- log(dgamma(1 / lb.inv, 1, scale = gammab))
  
  return(l1 + l2 + l3)
  
}


