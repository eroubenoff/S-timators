#------------------------------------------------------------------------------#

##### Shift Estimators  #####
##### Ethan Roubenoff   #####
##### Test No. 1, 7 Jan 2020 #####

# This script is an implementation of Tom Cassidy's Shift estimator (S-Timator).
# Given two incomplete mortality curves for two independent populations (men 
# women, for ex), the S-timator estimates the _difference_ in life expectancy 
# between the two.  The main assumption is:
#         m_male(a) = m_female(a + S)
# that one group's mortality is a linear scaling of the other's.

# Note that this model is implemented in discrete time, and correct conversions
# between integrals and summations should be verified.
# Required external data sources: two files simulating men and women's
# single-year mortality.  PDF of description is included. 

# In this script I am seeking to replicate Tom's work.  Next steps include the
# possibility of a spline, or a bayesian application.  I will also run this
# estimator on as many curves as possible. 

# A note on rounding: one of the big procedures here is accessing values at 
# given ages (m(70) = .08); ages are discrete here.  However, since our 
# estimator is not always an integer, there are inconsistencies with 
# m(alpha + delta) => m(20 + 8.2); there is no entry for m(28.2) in the 
# discrete data.  As a result, I just round to the nearest integer (m(28)).
# This may introduce systemic bias!

# Note that if any of the estimators are outside of the bounds of data,
# the function will break.  I do not have a good solution for this as of yet. 
#------------------------------------------------------------------------------#

##### Load data #####
library(tidyverse)
library(magrittr)
setwd("~/S-timators")
rm(list=ls())
gc()
m.sim <- readxl::read_excel("men-sim-file.xls", col_names = F) %>% 
  rename(age = "...1", v = "...2") %>%
  mutate(n = 1)
f.sim <- readxl::read_excel("women-sim-file.xls", col_names = F) %>%
  rename(age = "...1", v = "...2") %>%
  mutate(n = 1)

# note that the simulated data are in 1-year intervals. I have added a "n"
# column that can be used for non-constant age intervals.  These three
# columns (age, v, n) are the only necessary data.

m.sim
f.sim

ggplot() + 
  geom_line(data = m.sim, aes(x = age, y = v), color = "blue") + 
  geom_line(data = f.sim, aes(x = age, y = v), color = "red")  

# It can be seen that the female curve is right-shifted and scaled
# down relative to the male curve


#------------------------------------------------------------------------------#

##### Implementation #####
# Notation may be slightly different from the pdf
# all functions have prefix f_{quantity}
# CONSTANTS are in CAPITALS
# quantities should be consistent like: {variable}.{superscript}_{sex} sigma.sq.hat_m

mu_m <- weighted.mean(m.sim$age, m.sim$v)
mu_f <- weighted.mean(f.sim$age, f.sim$v)

ALPHA <- min(min(m.sim$age), min(f.sim$age))
BETA <- max(max(m.sim$age), max(f.sim$age))
DELTA <- round(mu_f - mu_m)

#### Iteration ####
# This function is a wrapper that returns all the constants needed for the
# estimators.  The default value of DELTA is mu_f - mu_m. 

stimator <- function(DELTA, ALPHA, BETA, f.sim, m.sim){
  
  if(ALPHA + DELTA < min (m.sim$age)) stop("ALPHA + DELTA falls outside of range")
  if(BETA - DELTA > max(m.sim$age)) stop("BETA - DELTA falls outside of range")
  
  #### mu ####
  browser()
  f_mu.hat <- function(alpha, beta, delta, df) {
    df %>% 
      filter(age > alpha & age < beta) %>%
      mutate(av = age * v * n ) %>%
      select(av) %>%
      sum() %>% 
      # divide_by(df %>% select(n) %>% sum()) %>%
      return()
  }
  
  # Note: need to round bc of discrete time
  mu.hat_w <- f_mu.hat(alpha = round(ALPHA+DELTA), beta = BETA, delta = DELTA, df = f.sim)
  mu.hat_m <- f_mu.hat(alpha = ALPHA, beta = round(BETA-DELTA), delta = DELTA, df = m.sim)
  
  # Quick sanity check:
  # ggplot() + 
  #   geom_line(data = m.sim, aes(x = age, y = v), color = "blue") + 
  #   geom_line(data = f.sim, aes(x = age, y = v), color = "red")  +
  #   geom_vline(xintercept = ALPHA + DELTA) + 
  #   geom_vline(xintercept = BETA - DELTA) + 
  #   geom_vline(xintercept = mu_m, color = "blue") + 
  #   geom_vline(xintercept = mu_f, color = "red") 
  
  #### sigma.hat.sq ####
  f_sigma.hat.sq <- function(alpha, beta,df, mu.hat) {
    df %>%
      filter(age > alpha, age < beta) %>%
      #mutate(a.mu = age - mu.hat) %>% # CHECK THIS
      mutate(a.mu = age) %>%
      mutate(r = a.mu^2 * v * n) %>%
      select(r) %>%
      sum() %>%
      return()
  }
  
  sigma.hat.sq_w <- f_sigma.hat.sq(alpha = round(ALPHA + DELTA), 
                                   beta = BETA, 
                                   df = f.sim, 
                                   mu.hat = mu.hat_w)
  
  sigma.hat.sq_m <- f_sigma.hat.sq(alpha = ALPHA, 
                                   beta = round(BETA - DELTA), 
                                   df = m.sim, 
                                   mu.hat = mu.hat_m)
  
  #### Q and Q.tilde ####
  f_Q <- function(alpha, beta, df) {
    df %>%
      filter(age > alpha, age < beta) %>%
      mutate(r = v*n) %>%
      sum()
  }
  
  Q <- f_Q(alpha = round(ALPHA + DELTA), 
           beta = BETA, 
           df = f.sim)
  Q.tilde <- f_Q(alpha = ALPHA, 
           beta = round(BETA - DELTA), 
           df = m.sim)
  
  #### R and R.tilde ####
  
  R <- mu.hat_m + DELTA
  R.tilde <- mu.hat_w - DELTA
  
  
  #### N and N.tilde ####
  f_N <- function(alpha, beta, R, sigma.hat.sq_m, Q, 
                  Q.tilde, df, mu.hat_m, mu.hat_w, delta ){
    a1 <- ((beta - R)^2 - sigma.hat.sq_m/Q.tilde) * (df %>% filter(age == beta) %>% pull("v"))
    a2 <- ((alpha - mu.hat_m)^2 - sigma.hat.sq_m/Q.tilde) * (df %>% filter(age == round(alpha + delta)) %>% pull("v"))
    a3 <- -2*mu.hat_w + 2*R*Q
    return(a1-a2+a3)
    }
  N <- f_N(ALPHA, BETA, R, sigma.hat.sq_m, Q,
           Q.tilde, f.sim,mu.hat_m, mu.hat_w, DELTA)  
  
  f_N.tilde <- function(alpha, beta, R.tilde, sigma.hat.sq_w, 
                        Q, Q.tilde, df, mu.hat_m, mu.hat_w, delta ){
    a1 <- ((beta - mu.hat_w)^2 - sigma.hat.sq_w/Q) *(df %>% filter(age == round(beta - delta)) %>% pull("v"))
    a2 <- ((alpha - R.tilde)^2 - sigma.hat.sq_w/Q)* (df %>% filter(age == alpha) %>% pull("v"))
    a3 <- -2*mu.hat_m + 2*R.tilde*Q.tilde
    return(a1 - a2 + a3)
  }
  
  N.tilde <- f_N.tilde(ALPHA, BETA, R.tilde, sigma.hat.sq_w, 
                       Q, Q.tilde, m.sim, mu.hat_m, mu.hat_w, DELTA)
  
  
  #### D and D.tilde ####
  # Need to calc v' (vp; v prime)
  # update this to use diff() function
  f.sim$vp <- 0
  m.sim$vp <- 0
  for (i in 1:nrow(f.sim)){
    vp <- (f.sim[i+1, "v"] - f.sim[i, "v"] ) %>% pull()
    n <- (f.sim[i, "n"]) %>% pull
    f.sim[i, "vp"] <- vp/n
  }
  f.sim[nrow(f.sim), "vp"] <- 0
  
  for (i in 1:nrow(m.sim)){
    vp <- (m.sim[i+1, "v"] - m.sim[i, "v"] ) %>% pull()
    n <- (m.sim[i, "n"]) %>% pull
    m.sim[i, "vp"] <- vp/n
  }
  m.sim[nrow(m.sim), "vp"] <- 0
  
  f_D <- function(beta, R, sigma.hat.sq_m, Q.tilde, 
                  df, alpha, mu.hat_m, delta, Q){
    a1 <- ((beta - R)^2 - sigma.hat.sq_m/Q.tilde) * (df %>% filter(age == beta) %>% pull("vp"))
    a2 <- ((alpha - mu.hat_m)^2 - sigma.hat.sq_m/Q.tilde) * (df %>% filter(age == round(alpha + delta)) %>% pull("vp"))
    a3 <- 2*(beta - R)*(df %>% filter(age == beta) %>% pull("v")) 
    a4 <- 2*(alpha - mu.hat_m)*(df %>% filter(age == round(alpha + delta)) %>% pull("v")) + 2*Q
    return(a1 - a2 - a3 + a4)
  }
  
  D <- f_D(BETA, R, sigma.hat.sq_m, Q.tilde, m.sim, ALPHA, mu.hat_m, DELTA, Q)
  
  f_D.tilde <- function(beta, mu.hat_w, sigma.hat.sq_w, Q, df, alpha, R.tilde, Q.tilde, delta) {
    a1 <- ((beta - mu.hat_w)^2 - sigma.hat.sq_w/Q) * (df %>% filter(age == round(beta - delta)) %>% pull("vp"))
    a2 <- ((alpha - R.tilde)^2 - sigma.hat.sq_w/Q) * (df %>% filter(age == alpha) %>% pull("vp"))
    a3 <- 2 * (beta - mu.hat_w) * (df %>% filter(age == round(beta - delta)) %>% pull("v"))
    a4 <- 2*(alpha - R.tilde) * (df %>% filter(age == alpha) %>% pull("v")) + 2*Q.tilde
    return(a1 - a2 - a3 + a4)
  }
  
  D.tilde <- f_D.tilde(BETA, mu.hat_w, sigma.hat.sq_w, Q, m.sim, ALPHA, R.tilde, Q.tilde, DELTA)
  
  
  ret <- list("ALPHA" = ALPHA,
              "BETA" = BETA,
              "DELTA" = DELTA,
              "D" = D,
              "D.tilde" = D.tilde,
              "mu_f" = mu_f,
              "mu_m" = mu_m,
              "mu.hat_m" = mu.hat_m,
              "mu.hat_w" = mu.hat_w,
              "N" = N,
              "N.tilde" = N.tilde,
              "Q" = Q,
              "Q.tilde" = Q.tilde,
              "R" = R,
              "R.tilde" = R.tilde,
              "sigma.hat.sq_m" = sigma.hat.sq_m,
              "sigma.hat.sq_w" = sigma.hat.sq_w)
  
  return(ret)
}

ret <- stimator(DELTA, ALPHA, BETA, f.sim, m.sim)

#### The Formulae  ####
#### Estimator 1 ####
f_e1 <- function(delta = DELTA, 
                 N = N, 
                 D = D, 
                 sigma.hat.sq_w = sigma.hat.sq_w, 
                 sigma.hat.sq_m = sigma.hat.sq_m, 
                 mu.hat_w = mu.hat_w, 
                 Q = Q, 
                 Q.tilde = Q.tilde, 
                 R = R) {
  a1 <- delta - N/D
  a2 <- N^2 - 2*D*(sigma.hat.sq_w + 2*mu.hat_w^2 - Q*mu.hat_w^2 - 2*R*mu.hat_w + 
                     Q*R^2 -  Q/Q.tilde*sigma.hat.sq_m)
  browser()
  return(a1 - (sqrt(a2)/D))
}

e1 <- with(ret, f_e1(delta = DELTA, 
           N = N,
           D = D, 
           sigma.hat.sq_w = sigma.hat.sq_w, 
           sigma.hat.sq_m = sigma.hat.sq_m, 
           mu.hat_w = mu.hat_w, 
           Q = Q, 
           Q.tilde = Q.tilde, 
           R = R) )

print(e1)


#### Estimator 2 ####
f_e2 <- function(delta = DELTA, 
                 N.tilde = N.tilde,
                 D.tilde = D.tilde,
                 sigma.hat.sq_m = sigma.hat.sq_m,
                 mu.hat_m = mu.hat_m,
                 Q.tilde = Q.tilde,
                 R.tilde = R.tilde,
                 Q = Q,
                 sigma.hat.sq_w = sigma.hat.sq_w) {
  
  a1 <- delta + N.tilde/D.tilde
  a2 <-  N.tilde^2 - 2*D.tilde * (sigma.hat.sq_m + 2*mu.hat_m^2 - 
                                  Q.tilde*mu.hat_m^2 - 2*R.tilde*mu.hat_m +
                                    Q.tilde * R.tilde^2 - Q.tilde/Q*sigma.hat.sq_w)
  
  return(a1 + sqrt(a2)/D.tilde)
  
}
e2 <- with(ret, f_e2(delta = DELTA, 
           N.tilde = N.tilde,
           D.tilde = D.tilde,
           sigma.hat.sq_m = sigma.hat.sq_m,
           mu.hat_m = mu.hat_m,
           Q.tilde = Q.tilde,
           R.tilde = R.tilde,
           Q = Q,
           sigma.hat.sq_w = sigma.hat.sq_w))
print(e2)

#### Estimator 3 ####

f_e3 <- function(delta = DELTA,
                 sigma.hat.sq_w = sigma.hat.sq_w,
                 mu.hat_w = mu.hat_w,
                 Q = Q,
                 R = R, 
                 Q.tilde = Q.tilde,
                 sigma.hat.sq_m = sigma.hat.sq_m,
                 N=N) {
  
  a1 <- sigma.hat.sq_w + 2*mu.hat_w^2 - Q.tilde*mu.hat_w^2 - 2*R*mu.hat_w +
    Q.tilde*R^2 - Q/Q.tilde*sigma.hat.sq_m
  return(delta - a1/N)
}
 
e3 <- with(ret, f_e3(delta = DELTA,
           sigma.hat.sq_w = sigma.hat.sq_w,
           mu.hat_w = mu.hat_w,
           Q = Q,
           R = R, 
           Q.tilde = Q.tilde,
           sigma.hat.sq_m = sigma.hat.sq_m,
           N=N))
print(e3)

#### Estimator 4 ####

f_e4 <- function(delta = DELTA,
                 sigma.hat.sq_m = sigma.hat.sq_m,
                 mu.hat_m= mu.hat_m,
                 Q.tilde = Q.tilde,
                 R.tilde = R.tilde,
                 Q = Q,
                 sigma.hat.sq_w = sigma.hat.sq_w,
                 N.tilde = N.tilde) {
  a1 <- sigma.hat.sq_m + 2*mu.hat_m^2 - Q.tilde*mu.hat_m^2 - 2*R.tilde*mu.hat_m +
    Q.tilde * R.tilde^2 - Q.tilde/Q*sigma.hat.sq_w
  return(delta + a1/N.tilde)
}

e4 <- with(ret, f_e4(delta = DELTA,
           sigma.hat.sq_m = sigma.hat.sq_m,
           mu.hat_m= mu.hat_m,
           Q.tilde = Q.tilde,
           R.tilde = R.tilde,
           Q = Q,
           sigma.hat.sq_w = sigma.hat.sq_w,
           N.tilde = N.tilde))
  
print(e4)

#### Iteration Procedure ####
# Program has already run once, providing first estimates of e1, e2, e3, e4
# Create nx4 matrix of values
# Call iteratively replacing DELTA with e1
# After completing n iterations, run to clear, and repeat with e2, e3, e4

n <- 10
iteration_df <- data.frame(matrix(0, n, 4))
colnames(iteration_df) <- c("e1", "e2", "e3", "e4")

# e1
iteration_df[1, "e1"] <- e1
for (i in 2:n) {
  ret <- stimator(DELTA = e1, ALPHA, BETA, f.sim, m.sim)
  e1 <- with(ret, f_e1(delta = e1,
                       N = N, 
                       D = D, 
                       sigma.hat.sq_w = sigma.hat.sq_w, 
                       sigma.hat.sq_m = sigma.hat.sq_m, 
                       mu.hat_w = mu.hat_w, 
                       Q = Q, 
                       Q.tilde = Q.tilde, 
                       R = R))
  iteration_df[i, "e1"] <- e1
}
iteration_df
plot(iteration_df$e3)

# e4
iteration_df[1, "e4"] <- e4
for (i in 2:n) {
  ret <- stimator(DELTA = e4, ALPHA, BETA, f.sim, m.sim)
  e4 <- with(ret, f_e4(delta = DELTA,
             sigma.hat.sq_m = sigma.hat.sq_m,
             mu.hat_m= mu.hat_m,
             Q.tilde = Q.tilde,
             R.tilde = R.tilde,
             Q = Q,
             sigma.hat.sq_w = sigma.hat.sq_w,
             N.tilde = N.tilde))
  iteration_df[i, "e4"] <- e4
}

# e3
iteration_df[1, "e3"] <- e3
for (i in 2:n) {
  ret <- stimator(DELTA = e3, ALPHA, BETA, f.sim, m.sim)
  e3 <- with(ret, f_e3(delta = e3,
                       sigma.hat.sq_w = sigma.hat.sq_w,
                       mu.hat_w = mu.hat_w,
                       Q = Q,
                       R = R, 
                       Q.tilde = Q.tilde,
                       sigma.hat.sq_m = sigma.hat.sq_m,
                       N=N))
  iteration_df[i, "e3"] <- e3
}
iteration_df
plot(iteration_df$e3)

# e4
iteration_df[1, "e4"] <- e4
for (i in 2:n) {
  ret <- stimator(DELTA = e4, ALPHA, BETA, f.sim, m.sim)
  e4 <- with(ret, f_e4(delta = DELTA,
             sigma.hat.sq_m = sigma.hat.sq_m,
             mu.hat_m= mu.hat_m,
             Q.tilde = Q.tilde,
             R.tilde = R.tilde,
             Q = Q,
             sigma.hat.sq_w = sigma.hat.sq_w,
             N.tilde = N.tilde))
  iteration_df[i, "e4"] <- e4
}
iteration_df
plot(iteration_df$e4)

















  