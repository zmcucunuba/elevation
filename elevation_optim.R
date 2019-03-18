rm(list=ls())

library(readxl)
dat <- read_excel("Illustration for Kelly - estimating the f function.xlsx", 'Sheet1')
dat <- dat[1:1120,]
dat$trials <- 1

mod1 <- function(params)
{
  par_a  <- params[1]
  par_b  <- params[2]
  par_c  <- params[3]
  prob_zero <- params[4]
  
  z <- as.matrix(dat$elevation)
  pos <- as.matrix(dat$infected)
  trials <- as.matrix(rep(1, length(pos)))
  
  # prob <- ((1+exp(par_a))/(1+exp(par_a*(1-z/z_half))))*(1-par_c) + par_c
  prob <-  (( 1+exp( -par_a* par_b ) ) / (1+exp( par_a *( z - par_b ) ))  * (1-par_c) + par_c) * prob_zero
  
  LL <- (dbinom(x=pos, size=trials, prob=prob, log = TRUE))
  
  
  
  sum_LL <- -sum(LL)
  
  print(sum_LL)
  
  return(sum_LL)
  
} 


# params <- c(0.004858208, 1688.698703, 0.02615773, 0.491040964)
params_init <- c(0.1, 1600, 0.06, 0.4)

mod1(params_init)




res_optim <- optim(par = params_init,
                   fn = mod1,
                   NULL,
                   method = "L-BFGS-B",
                   lower=c(0, 1000, 0.001, 0.0001),
                   upper = c(1, 2000, 1, 0.9999))

round(res_optim$par, 4)


res_optim2 <- optim(par = res_optim$par,
                    fn = mod1,
                    NULL,
                    method = "L-BFGS-B",
                    lower=c(0, 1000, 0.001, 0.0001),
                    upper = c(1, 2000, 1, 0.9999))

round(res_optim2$par, 4)

res_optim2<- optim(par = params_model,
                   fn = mod1,
                   NULL,
                   method = "L-BFGS-B",
                   lower=c(0, 1000, 0.001, 0.0001),
                   upper = c(1, 2000, 1, 0.9999))


res_optim2$par
