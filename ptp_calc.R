##################################################################################################
#### CODES TO CALCULATE POSTERIOR TAIL PROBABILITY FOR THE SIMULATION STUDY IN THE PAPER
##################################################################################################

## SETTINGS 1 and 3

# homdir <- "H:/lenovo backup/staff scientist/combination prevention/codes/backup from loaner/cluster/workspaces/batch_10_10_16_reimagined"
# setwd(homdir)

library(coda)
library(reshape2)


### CALCULATING PTP FOR SETTINGS 1 AND 3

MCMCsamples = 5000;
lagsize = 4
thetasize = 17
nloc = 6

index.vector <- function(sample, trace) seq(1, sample, trace)

ind_2 <- index.vector(MCMCsamples, 2)
ind_5 <- index.vector(MCMCsamples, 5)
ind_10 <- index.vector(MCMCsamples, 10)

bb = c('000000', '151515151515')
cc = c('TRUE', 'FALSE')

n.script = 20

B <- 100
K <- 5000 / B
meanP1 <- array(0, c(length(B), 2, 2))
for (j in 1:length(B)){
  P_MD1 <- array(NA, c(2, 2, n.script))
  for (k in 1:length(bb)){
    for (kk in 1:length(cc)){
      for (t in 1:n.script){
        tryCatch({load.nm  <- paste("MCMC_b_", bb[k], "_ut_", cc[kk], "_s_", t, ".Rdata", sep = "")
        load(load.nm)
        datmat <- mcmc.run$mcmc.chain[, 11:16]
        
        u_hat <- colMeans(datmat)
        uk_hat <- matrix(0, K[j], length(u_hat))
        m_hat <- array(0, c(K[j], length(u_hat), length(u_hat)))
        for (i in 1:K[j]){
          uk_hat[i, ] <- colMeans(datmat[((i-1)*B[j]+1):(i*B[j]), ])
          m_hat[i, , ] <- (uk_hat[i, ]-u_hat)%*%t(uk_hat[i, ]-u_hat)
        }
        tau_sq <- B[j]/(K[j]-1)*apply(m_hat, c(2, 3), sum)
        var.u <- tau_sq
        
        MDdat <- mahalanobis(datmat, u_hat, var.u)
        MDmean <- mahalanobis(u_hat, rep(0, 6), var.u)
        P_MD1[k, kk, t] <- mean(MDdat >= MDmean)
        }, error = function(e){print('no such WS')})
      }
    }
  }
  meanP1[j, , ] <- apply(P_MD1, c(1, 2), function(x) mean(x, na.rm = T))
}

print(meanP1)


### CALCULATING PTP FOR SETTINGS 2

# homdir <- "H:/lenovo backup/staff scientist/combination prevention/codes/backup from loaner/new_codes/batch_5_31_17/onetheta_extra_04_30_18"
# setwd(homdir)
MCMCsamples = 4000;

bb = c('flat', 'narrow')
cc = c('TRUE', 'FALSE')

B <- 100
K <- 4000 / B
meanP1 <- array(0, c(length(B), 2, 2))
for (j in 1:length(B)){
  P_MD1 <- array(NA, c(2, 2, n.script))
  for (k in 1:length(bb)){
    for (kk in 1:length(cc)){
      for (t in 1:n.script){
        tryCatch({load.nm  <- paste("onetheta_", bb[k], "_ind_5_tb_-50_ut_", cc[kk], "_s_", t, ".Rdata", sep = "")
        load(load.nm)
        datmat <- mcmc.run$mcmc.chain[, 2:7]
        
        u_hat <- colMeans(datmat)
        uk_hat <- matrix(0, K[j], length(u_hat))
        m_hat <- array(0, c(K[j], length(u_hat), length(u_hat)))
        for (i in 1:K[j]){
          uk_hat[i, ] <- colMeans(datmat[((i-1)*B[j]+1):(i*B[j]), ])
          m_hat[i, , ] <- (uk_hat[i, ]-u_hat)%*%t(uk_hat[i, ]-u_hat)
        }
        tau_sq <- B[j]/(K[j]-1)*apply(m_hat, c(2, 3), sum)
        var.u <- tau_sq
        
        MDdat <- mahalanobis(datmat, u_hat, var.u)
        MDmean <- mahalanobis(u_hat, rep(0, 6), var.u)
        P_MD1[k, kk, t] <- mean(MDdat >= MDmean)
        }, error = function(e){print('no such WS')})
      }
    }
  }
  meanP1[j, , ] <- apply(P_MD1, c(1, 2), function(x) mean(x, na.rm = T))
}

print(meanP1)
