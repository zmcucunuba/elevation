
rm(list=ls())

library(readxl)
dat <- read_excel("Illustration for Kelly - estimating the f function.xlsx", 'Sheet1')
dat <- dat[1:1120,]
dat$trials <- 1

mod2 <- function(params)
{
  par_a  <- params[1]
  z_half <- params[2]
  par_c  <- params[3]
  
  z <- as.matrix(dat$elevation)
  pos <- as.matrix(dat$infected)
  trials <- as.matrix(rep(1, length(pos)))
  
  prob <- ((1+exp(par_a))/(1+exp(par_a*(1-z/z_half))))*(1-par_c) + par_c
  
  LL <- (dbinom(x=pos, size=trials, prob=prob, log = TRUE))
  
  if(any(is.nan(LL)))
     {
    
    browser()
    
  }
  
  sum_LL <- sum(LL)
  
  return(sum_LL)
  
} 

params <- c(-10.1, 1703.935949, 0.03594105)

mod2(params)

res_optim <- optim(par = params,
      fn = mod2,
      NULL,
      method = "L-BFGS-B",
      lower=c(-20, 1000, 0.001),
      upper = c(1, 3500, 0.1))


if(res_optim$convergence != 0)
{
  res_optim <- optim(par = params,
                     fn = mod2,
                     NULL,
                     method = "L-BFGS-B",
                     lower=c(-20, 1000, 0.001),
                     upper = c(1, 3500, 0.1))
  
}
res_optim$convergence 


prob <- function(params)
{
  par_a  <- params[1]
  z_half <- params[2]
  par_c  <- params[3]
  
  z <- as.matrix(dat$elevation)
  pos <- as.matrix(dat$infected)
  trials <- as.matrix(rep(1, length(pos)))
  
  prob <- ((1+exp(par_a))/(1+exp(par_a*(1-z/z_half))))*(1-par_c) + par_c
  
  return(prob)
  
} 

params_christl <- c(-10.1,
                  1703.9359492813,
                  0.0359410497799828)

params_optim <- res_optim$par
res <- prob(params_optim)
dat$prob_model <- res
plot(dat$elevation, dat$infected)
lines(dat$elevation, dat$prob_model, col= 'red')

