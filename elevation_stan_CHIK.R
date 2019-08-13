rm(list=ls())
library(rstan)
library(ggplot2) #Opening libraries
library(rstan)
library(reshape2)
library(dplyr)
library(bayesplot)
library(loo)
library(rstanarm)
library(rstan)
library(readxl)
library(grid)
library(gridExtra)
library(readr)
library(tidyverse)

# elevation_model_stan  <-  stan_model('elevation_model.stan')
# saveRDS(elevation_model_stan, "elevation_model_stan.RDS")
elevation_model_stan  <-  readRDS('elevation_model_stan.RDS')


dat  <- readRDS('data/elevation_infected_chik.RDS')

dat$trials <- 1

stan_data <- list(
  Nobs    = nrow(dat),
  Npos    = dat$infected,
  Ntrials = dat$trials,
  max_elevation = max(dat$elevation),
  z       = dat$elevation)


#  Run the Model
n_iter <- 8000
fit <- sampling(elevation_model_stan, 
                data=stan_data, 
                iter = n_iter)

parameters <- c("a", "b", 'c', 'p_zero')
burnin <- floor(n_iter/2)
res_total <- data.frame(fit@sim$samples)[(burnin+1):n_iter,]
res_total <- res_total [,parameters]


plot_info_table <- function(info){
  blank <- data.frame(x= 1:10, y = 1:10)
  pi <- ggplot(blank, aes(x, y)) +
    geom_blank() + ylab('') + xlab ('') + theme_void(10) +
    annotation_custom(tableGrob(info))
  return(pi)
}




pdf('res_chik_corrected.pdf')
plot(fit, plotfun = "trace", pars = parameters, inc_warmup = TRUE)
par(mfrow=c(2,2))
hist(res_total$a, col = 'red', main = 'a')
hist(res_total$b, col = 'red',  main = 'b')
hist(res_total$c, col = 'red',  main = 'c')
hist(res_total$p_zero, col = 'red',  main = 'p_zero')

results <- data.frame(parameter = parameters, lower = NA, median = NA, upper = NA)
results[1,2:4] <- data.frame(t(quantile(res_total$a, c(0.025, 0.5, 0.975))))
results[2,2:4] <- data.frame(t(quantile(res_total$b, c(0.025, 0.5, 0.975))))
results[3,2:4] <- data.frame(t(quantile(res_total$c, c(0.025, 0.5, 0.975))))
results[4,2:4] <- data.frame(t(quantile(res_total$p_zero, c(0.025, 0.5, 0.975))))

plot_info_table(results) + 
  ggtitle('Median and CrI of posterior dist')


write.csv(results, 'results_final_chik_corrected.csv')

estimated_pars <- results[,'median']


mod1 <- function(params)
{
  par_a  <- params[1]
  par_b  <- params[2]
  par_c  <- params[3]
  prob_zero <- params[4]
  
  z <- as.matrix(dat$elevation)
  pos <- as.matrix(dat$infected)
  trials <- as.matrix(rep(1, length(pos)))
  
  prob <-  (( 1+exp( -par_a* par_b ) ) / (1+exp( par_a *( z - par_b ) ))  * (1-par_c) + par_c) * prob_zero
  LL <- (dbinom(x=pos, size=trials, prob=prob, log = TRUE))
  sum_LL <- sum(LL)
  
  return(sum_LL)
  
} 




res_total$log_lik <- NA
for ( i in 1:nrow(res_total)) {
  test_pars <- c(res_total$a[i],
                 res_total$b[i],
                 res_total$c[i],
                 res_total$p_zero[i])
  res_total$log_lik[i] <- mod1(test_pars)
  
}

max_LL <- max(res_total$log_lik)
pars_max_LL <- res_total[res_total$log_lik == max_LL,]

pars_max_LL <- t(pars_max_LL)
colnames(pars_max_LL) <- 'Values at max LL'
plot_info_table(pars_max_LL) 


estimated_pars <- read_csv("results_final_chik_corrected.csv")

prob_mod <- function(params)
{
  par_a  <- params[1]
  par_b  <- params[2]
  par_c  <- params[3]
  prob_zero <- params[4]
  z <- as.matrix(dat$elevation)
  
  
  prob <-  (( 1+exp( -par_a* par_b ) ) / (1+exp( par_a *( z - par_b ) ))  * (1-par_c) + par_c) 
  return(prob)
  
} 


par(mfrow=c(1,1))

dat$post_median <- prob_mod(estimated_pars$median)
dat$post_lower <- prob_mod(estimated_pars$lower)
dat$post_upper <- prob_mod(estimated_pars$upper)

dat <- arrange(dat, elevation)

plot(dat$elevation, dat$post_median, ylim = c(0,1), type = 'l', cex = 2)
lines(dat$elevation, dat$post_lower, col = 'blue')
lines(dat$elevation, dat$post_upper, col = 'blue')


dev.off()
