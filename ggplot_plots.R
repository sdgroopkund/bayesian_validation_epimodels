######################################################
#### CODES TO GENERATE THE PLOTS IN THE PAPER
######################################################

library(reshape2)
library(ggplot2)
library(ggdendro)
library(cowplot)
library(dotwhisker)
library(dplyr)
library(broom)
library(ggsignif)
library(contrast)
library(coda)

########################################################################################################################################################

## DISCREPANCY PLOTS (SETTINGS 1 and 3)
 
# homdir <- "H:/lenovo backup/staff scientist/combination prevention/codes/backup from loaner/cluster/workspaces/batch_10_10_16_reimagined"
# setwd(homdir)

## null scenario

bb = c('000000', '151515151515')
cc = c('TRUE', 'FALSE')

MCMCsamples = 5000;
lagsize = 4
thetasize = 17
nloc = 6

index.vector <- function(sample, trace) seq(1, sample, trace)

ind_2 <- index.vector(MCMCsamples, 2)
ind_5 <- index.vector(MCMCsamples, 5)
ind_10 <- index.vector(MCMCsamples, 10)

n.script = 20

plot1 <- list()
for (i in 1:nloc){
  p <- NULL
  for (t in 1:n.script){
    tryCatch({load.nm  <-  paste("MCMC_b_", bb[1], "_ut_", cc[2], "_s_", t, ".Rdata", sep = "")    
    load(load.nm)
    mcmc.obj1 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data1 <- data.frame('bias' = as.numeric(mcmc.obj1[, (10 + i)]))}, error = function(e){print('no such WS')})
    
    tryCatch({load.nm  <- paste("MCMC_b_", bb[1], "_ut_", cc[1], "_s_", zz[t], ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj2 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data2 <- data.frame('bias' = as.numeric(mcmc.obj2[, (10 + i)]))}, error = function(e){print('no such WS')})
    
    if (t == 1) {p <- ggplot(mcmc.data1,  aes(x = bias))  +  
      stat_density(size = 0.75, aes(x = bias, color = 'brown'), geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2, aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity') + 
      theme_minimal() + coord_cartesian(xlim  =  c(-0.5,  0.5),  ylim  =  c(0, 9)) 
    }
    else if (t > 1) {p <- p + stat_density(data = mcmc.data1,  aes(x = bias,  color = 'brown'), size = 0.75, geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2,  aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity')}
  }
  q <- p + ggtitle(paste('density at location ', i, sep = "")) + xlab("") + ylab("") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)) + 
    theme(axis.text = element_text(size = 15, face = 'bold')) + 
    theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
    theme(panel.background  =  element_blank()) + 
    scale_y_continuous(expand  =  c(0,  0)) + 
    scale_colour_manual(name  =  'Update Prior',  
                        values  = c('brown' = 'brown', 'steelblue' = 'steelblue'),  labels  =  c('No', 'Yes')) + 
    theme(legend.position  =  c(0.8, 0.8) ) + guides(color  =  guide_legend(override.aes  =  list(size = 1.5)))
  
  
  plot1[[i]] <- q
}

dev.new(width = 1700,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot1 <- plot_grid(plot1[[1]], plot1[[2]], plot1[[3]], plot1[[4]], plot1[[5]], plot1[[6]],  ncol = 3)

title  <-  ggdraw()  +  draw_label("posterior distribution of discrepancy under the null scenario",  fontface = 'bold', size = 18,  hjust = 0.5)

FP1.new <- plot_grid(title,  bigplot1,  ncol = 1,  rel_heights = c(0.1,  1))
FP1.new



## faulty structure

plot2 <- list()
for (i in 1:nloc){
  p <- NULL
  for (t in 1:n.script){
    tryCatch({load.nm  <-  paste("MCMC_b_", bb[2], "_ut_", cc[2], "_s_", t, ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj1 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data1 <- data.frame('bias' = as.numeric(mcmc.obj1[, (10 + i)]))}, error = function(e){print('no such WS')})
    
    tryCatch({load.nm  <- paste("MCMC_b_", bb[2], "_ut_", cc[1], "_s_", t, ".Rdata", sep = "")    
    load(load.nm)
    mcmc.obj2 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data2 <- data.frame('bias' = as.numeric(mcmc.obj2[, (10 + i)]))}, error = function(e){print('no such WS')})
    
    if (t == 1) {p <- ggplot(mcmc.data1,  aes(x = bias))  +  
      stat_density(size = 0.75, aes(x = bias, color = 'brown'), geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2, aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity') + 
      theme_minimal() + coord_cartesian(xlim  =  c(-0.5,  0.5),  ylim  =  c(0, 9)) 
    }
    else if (t > 1) {p <- p + stat_density(data = mcmc.data1,  aes(x = bias,  color = 'brown'), size = 0.75, geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2,  aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity')}
  }
  q <- p + ggtitle(paste('density at location ', i, sep = "")) + xlab("") + ylab("") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)) + 
    theme(axis.text = element_text(size = 15, face = 'bold')) + 
    theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
    theme(panel.background  =  element_blank()) + 
    scale_y_continuous(expand  =  c(0,  0)) + 
    scale_colour_manual(name  =  'Update Prior',  
                        values  = c('brown' = 'brown', 'steelblue' = 'steelblue'),  labels  =  c('No', 'Yes')) + 
    theme(legend.position  =  c(0.87, 0.85) ) + guides(color  =  guide_legend(override.aes  =  list(size = 1.5)))
  
  
  plot2[[i]] <- q
}

dev.new(width = 1700,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot2 <- plot_grid(plot2[[1]], plot2[[2]], plot2[[3]], plot2[[4]], plot2[[5]], plot2[[6]],  ncol = 3)

title  <-  ggdraw()  +  draw_label("posterior distribution of discrepancy under faulty structure",  fontface = 'bold', size = 18,  hjust = 0.5)

FP2.new <- plot_grid(title,  bigplot2,  ncol = 1,  rel_heights = c(0.1,  1))
FP2.new

#########################################################################################################################################################

## DISCREPANCY PLOTS FOR SETTING 2

# homdir <- "H:/lenovo backup/staff scientist/combination prevention/codes/backup from loaner/new_codes/batch_5_31_17/onetheta_extra_04_30_18"
# setwd(homdir)

MCMCsamples = 4000;
lagsize = 4
thetasize = 17
nloc = 6
thetaind <- c(5)
n.script = 10
index.vector <- function(sample, trace) seq(1, sample, trace)

ind_2 <- index.vector(MCMCsamples, 2)
ind_5 <- index.vector(MCMCsamples, 5)
ind_10 <- index.vector(MCMCsamples, 10)

bb = c('flat', 'narrow')
cc = c('TRUE', 'FALSE') 

## FLAT PRIOR 
plot1 <- list()
for (i in 1:nloc){
  p <- NULL
  for (t in 1:n.script){
    tryCatch({load.nm  <- paste("onetheta_", bb[1], "_ind_5_tb_-50_ut_", cc[2], "_s_", t, ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj1 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data1 <- data.frame('bias' = as.numeric(mcmc.obj1[, (1 + i)]))}, error = function(e){print('no such WS')})
    
    tryCatch({load.nm  <- paste("onetheta_", bb[1], "_ind_5_tb_-50_ut_", cc[1], "_s_", t, ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj2 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data2 <- data.frame('bias' = as.numeric(mcmc.obj2[, (1 + i)]))}, error = function(e){print('no such WS')})
    
    if (t == 1) {p <- ggplot(mcmc.data1,  aes(x = bias))  +  
      stat_density(size = 0.75, aes(x = bias, color = 'brown'), geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2, aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity') + 
      theme_minimal() + coord_cartesian(xlim  =  c(-0.6,  0.6),  ylim  =  c(0, 5.5)) 
    }
    else if (t > 1) {p <- p + stat_density(data = mcmc.data1,  aes(x = bias,  color = 'brown'), size = 0.75, geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2,  aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity')}
  }
  q <- p + ggtitle(paste('density at location ', i, sep = "")) + xlab("") + ylab("") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)) + 
    theme(axis.text = element_text(size = 15, face = 'bold')) + 
    theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
    theme(panel.background  =  element_blank()) + 
    scale_y_continuous(expand  =  c(0,  0)) + 
    scale_colour_manual(name  =  'Update Prior',  
                        values  = c('brown' = 'brown', 'steelblue' = 'steelblue'),  labels  =  c('No', 'Yes')) + 
    theme(legend.position  =  c(0.8, 0.8) ) + guides(color  =  guide_legend(override.aes  =  list(size = 1.5)))
  
  
  plot1[[i]] <- q
}

dev.new(width = 1700,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot1 <- plot_grid(plot1[[1]], plot1[[2]], plot1[[3]], plot1[[4]], plot1[[5]], plot1[[6]],  ncol = 3)

title  <-  ggdraw()  +  draw_label(expression(bold(paste('posterior distribution of discrepancy for negative bias in ', sigma[1],  ' when prior is flat'))),  fontface = 'bold', size = 18,  hjust = 0.5)

FP3.new <- plot_grid(title,  bigplot1,  ncol = 1,  rel_heights = c(0.1,  1))
FP3.new




## For narrow prior

plot2 <- list()
for (i in 1:nloc){
  p <- NULL
  for (t in c(1:8, 10, 11)){
    tryCatch({load.nm  <- paste("onetheta_", bb[2], "_ind_5_tb_-50_ut_", cc[2], "_s_", t, ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj1 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data1 <- data.frame('bias' = as.numeric(mcmc.obj1[, (1 + i)]))}, error = function(e){print('no such WS')})
    
    tryCatch({load.nm  <- paste("onetheta_", bb[2], "_ind_5_tb_-50_ut_", cc[1], "_s_", t, ".Rdata", sep = "")
    load(load.nm)
    mcmc.obj2 <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
    mcmc.data2 <- data.frame('bias' = as.numeric(mcmc.obj2[, (1 + i)]))}, error = function(e){print('no such WS')})
    
    if (t == 1) {p <- ggplot(mcmc.data1,  aes(x = bias))  +  
      stat_density(size = 0.75, aes(x = bias, color = 'brown'), geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2, aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity') + 
      theme_minimal() + coord_cartesian(xlim  =  c(-0.6,  0.6),  ylim  =  c(0, 6)) 
    }
    else if (t > 1) {p <- p + stat_density(data = mcmc.data1,  aes(x = bias,  color = 'brown'), size = 0.75, geom = 'line', position = 'identity') + 
      stat_density(data = mcmc.data2,  aes(x = bias, color = 'steelblue'), size = 0.75, geom = 'line', position = 'identity')}
  }
  q <- p + ggtitle(paste('density at location ', i, sep = "")) + xlab("") + ylab("") + 
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 18)) + 
    theme(axis.text = element_text(size = 15, face = 'bold')) + 
    theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
    theme(panel.background  =  element_blank()) + 
    scale_y_continuous(expand  =  c(0,  0)) + 
    scale_colour_manual(name  =  'Update Prior',  
                        values  = c('brown' = 'brown', 'steelblue' = 'steelblue'),  labels  =  c('No', 'Yes')) + 
    theme(legend.position  =  c(0.8, 0.8) ) + guides(color  =  guide_legend(override.aes  =  list(size = 1.5)))
  
  
  plot2[[i]] <- q
}

dev.new(width = 1700,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot2 <- plot_grid(plot2[[1]], plot2[[2]], plot2[[3]], plot2[[4]], plot2[[5]], plot2[[6]],  ncol = 3)

title  <-  ggdraw()  +  draw_label(expression(bold(paste('posterior distribution of discrepancy for negative bias in ', sigma[1],  ' when prior is narrow'))),  fontface = 'bold', size = 18,  hjust = 0.5)

FP4.new <- plot_grid(title,  bigplot2,  ncol = 1,  rel_heights = c(0.1,  1))
FP4.new

##################################################################################################################################################################################################

#### PLOT FOR POSTERIOR DISTRIBUTION OF THETA 5
## Flat Prior

plot.flat <- list(2)
for (i in 2:1){  
  which.index = 1
  which.xlim = c(0, 0.5)
  load.nm  <- paste("onetheta_", bb[1], "_ind_5_tb_-50_ut_", cc[i], "_s_1.Rdata", sep = "")
  load(load.nm)
  prior.samples.f <- sapply(1:N,  function(i) logitgauss.prior(u.m[theta.ind], prec = delta.prec))
  mcmc.psf <- mcmc(prior.samples.f)
  mcmc.data2 <- data.frame('density' = as.numeric(mcmc.psf))
  p <- ggplot(mcmc.data2,  aes(x = density)) + stat_density(geom = 'line', aes(colour = 'brown'), size = 0.75) + 
    geom_vline(aes(xintercept = u.start[theta.ind], colour = 'black'), size = 1.5) + 
    theme_minimal() + coord_cartesian(xlim  =  c(0.027, 0.5),  ylim  =  c(0, 15))
  for (t in c(1:10)){
    tryCatch({
      load.nm  <- paste("onetheta_", bb[1], "_ind_5_tb_-50_ut_", cc[i], "_s_", t, ".Rdata", sep = "")
      load(load.nm)
      mcmc.obj <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
      mcmc.data <- data.frame('density' = as.numeric(mcmc.obj[, which.index]))
      p <- p + stat_density(data = mcmc.data,  geom = 'line',  aes(x = density,  colour = 'steelblue'), size = 0.75)
    }, error = function(e){print('no such WS')
    })
  }
  p <- p + stat_density(data = mcmc.data2,  geom = 'line',  aes(x = density,  colour = 'brown'),  size = 0.75)
  
  if (i == 1) {
    q <- p + ggtitle(expression(bold(paste('Distribution for ', sigma[1], ' when we update its prior (flat)')))) + xlab("") + ylab("") + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15)) + 
      theme(axis.text = element_text(size = 15, face = 'bold')) + 
      theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
      theme(panel.background  =  element_blank())
  } else if (i == 2) {
    q <- p + ggtitle(expression(bold(paste('Distribution for ', sigma[1], ' when we do not update its prior (flat)')))) + xlab("") + ylab("") + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15)) + 
      theme(axis.text = element_text(size = 15, face = 'bold')) + 
      theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
      theme(panel.background  =  element_blank())
  }
  
  r <- q + scale_colour_manual(name  =  expression(paste('Distribution for ', sigma[1])), 
                           values  = c('brown' = "brown",  'black' = "black",  'steelblue' = "steelblue"),  
                           labels  =  c('true value', 'prior distribution', 'final distribution in \ndifferent MCMC runs')) + 
    theme(legend.position  =   c(0.8, 0.7),  legend.text = element_text(size = 16),  legend.title  =  element_text(size = 18))
  
  plot.flat[[i]] <- r
  
  
  # dev.off()
}

dev.new(width = 1100,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot5 <- plot_grid(plot.flat[[2]], plot.flat[[1]],  ncol = 1,  labels = c('(i)', '(ii)'))
title1  <-  ggdraw()  +  draw_label(expression(bold(paste('A. Flat prior for ', sigma[1]))),  fontface = 'bold', size = 22,  hjust = 0.5)
FP5 <- plot_grid(title1,  bigplot5,  ncol = 1,  rel_heights = c(0.05,  1))
FP5

dev.new(width = 1100,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
title1  <-  ggdraw()  +  draw_label(expression(bold(paste('A. Flat prior for ', sigma[1]))),  fontface = 'bold', size = 22,  hjust = 0.5)
FP5.n <- plot_grid(title1,  plot.flat[[1]],  ncol = 1,  rel_heights = c(0.05,  1))
FP5.n


## Narrow Prior

plot.narrow <- list(2)
for (i in 2:1){  
  which.index = 1
  which.xlim = c(0, 0.5)
  load.nm  <- paste("onetheta_", bb[2], "_ind_5_tb_-50_ut_", cc[i], "_s_1.Rdata", sep = "")
  load(load.nm)
  prior.samples.f <- sapply(1:N,  function(i) logitgauss.prior(u.m[theta.ind], prec = delta.prec))
  mcmc.psf <- mcmc(prior.samples.f)
  mcmc.data2 <- data.frame('density' = as.numeric(mcmc.psf))
  p <- ggplot(mcmc.data2,  aes(x = density)) + stat_density(geom = 'line', aes(colour = 'brown'), size = 0.75) + 
    geom_vline(aes(xintercept = u.start[theta.ind], colour = 'black'), size = 1.5) + 
    theme_minimal() + coord_cartesian(xlim  =  c(0.027, 0.5),  ylim  =  c(0, 45))
  for (t in c(1:10)){
    tryCatch({
      load.nm  <- paste("onetheta_", bb[2], "_ind_5_tb_-50_ut_", cc[i], "_s_", t, ".Rdata", sep = "")
      load(load.nm)
      mcmc.obj <- mcmc(mcmc.run$mcmc.chain[ind_5, ])
      mcmc.data <- data.frame('density' = as.numeric(mcmc.obj[, which.index]))
      p <- p + stat_density(data = mcmc.data,  geom = 'line',  aes(x = density,  colour = 'steelblue'), size = 0.75)
    }, error = function(e){print('no such WS')
    })
  }
  p <- p + stat_density(data = mcmc.data2,  geom = 'line',  aes(x = density,  colour = 'brown'),  size = 0.75)
  
  if (i == 1) {
    q <- p + ggtitle(expression(bold(paste('Distribution for ', sigma[1], ' when we update its prior (narrow)')))) + xlab("") + ylab("") + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15)) + 
      theme(axis.text = element_text(size = 15, face = 'bold')) + 
      theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
      theme(panel.background  =  element_blank())
  } else if (i == 2) {
    q <- p + ggtitle(expression(bold(paste('Distribution for ', sigma[1], ' when we do not update its prior (narrow)')))) + xlab("") + ylab("") + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 15)) + 
      theme(axis.text = element_text(size = 15, face = 'bold')) + 
      theme(axis.line  =  element_line(colour  =  "black",  size = 1)) + 
      theme(panel.background  =  element_blank())
  }
  
  r <- q + scale_colour_manual(name  =  expression(paste('Distribution for ', sigma[1])), 
                           values  = c('brown' = "brown", 'black' = "black",   'steelblue' = "steelblue"),  
                           labels  =  c('true value', 'prior distribution', 'final distribution in \ndifferent MCMC runs')) + 
    theme(legend.position  =   c(0.8, 0.7),  legend.text = element_text(size = 16),  legend.title  =  element_text(size = 18))
  
  plot.narrow[[i]] <- r
  
  
  # dev.off()
}


dev.new(width = 1100,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
bigplot6 <- plot_grid(plot.narrow[[2]], plot.narrow[[1]],  ncol = 1,  labels = c('(i)', '(ii)'))
title2  <-  ggdraw()  +  draw_label(expression(bold(paste('B. Narrow prior for ', sigma[1]))),  fontface = 'bold', size = 22,  hjust = 0.5)

FP6 <- plot_grid(title2,  bigplot6,  ncol = 1,  rel_heights = c(0.05,  1))
FP6

dev.new(width = 1100,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
title2  <-  ggdraw()  +  draw_label(expression(bold(paste('B. Narrow prior for ', sigma[1]))),  fontface = 'bold', size = 22,  hjust = 0.5)
FP6.n <- plot_grid(title2,  plot.narrow[[1]],  ncol = 1,  rel_heights = c(0.05,  1))
FP6.n


dev.new(width = 2600,  height = 1600, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
FP.all  <-  plot_grid(title1,  title2,  bigplot5,  bigplot6,  ncol = 2,  rel_heights = c(0.05,  1))
FP.all

dev.new(width = 2600,  height = 1200, noRStudioGD  =  T)
theme_set(theme_cowplot(font_size = 12)) # reduce default font size
FP.all.n  <-  plot_grid(title1,  title2,  plot.flat[[1]],  plot.narrow[[1]],  ncol = 2,  rel_heights = c(0.05,  1))
FP.all.n

