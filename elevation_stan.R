

rm(list=ls())
library(rstan)
library(ggplot2) #Opening libraries
library(rstan)
library(reshape2)
library(dplyr)
library(bayesplot)
library(ggthemr)
library(loo)
library(rstanarm)
library(rstan)
library(readxl)



dat <- read_excel("Illustration for Kelly - estimating the f function.xlsx", 'Sheet1')
dat <- dat[1:1120,]
dat$trials <- 1

stan_data <- list(
  Nobs    = nrow(dat),
  Npos    = dat$infected,
  Ntrials = dat$trials,
  max_elevation = max(dat$elevation),
  z       = dat$elevation)

model_stan  <-  stanc("model.stan")
sm = stan_model(stanc_ret = model_stan, verbose=FALSE)

#  Run the Model
set.seed(123)
system.time(fit <- sampling(sm, 
                            data=stan_data, 
                            iter = 100))


plot(fit, plotfun = "trace", pars = c("a", "z_half", 'c'), inc_warmup = TRUE)



