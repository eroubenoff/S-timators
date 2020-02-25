#------------------------------------------------------------------------------#

##### Shift Estimators  #####
##### Ethan Roubenoff   #####
##### Continous Implementation: 6 Feb 2020 #####

# This script is an implementation of Tom Cassidy's Shift estimator (S-Timator).
# Given two incomplete mortality curves for two independent populations (men 
# women, for ex), the S-timator estimates the _difference_ in life expectancy 
# between the two.  The main assumption is:
#         m_male(a) = m_female(a + S)
# that one group's mortality is a linear scaling of the other's.

# This version is implemented in continuous time using linear interpolation, 
# which differs from the first other version (shift_discrete)
# Required external data sources: two files simulating men and women's
# single-year mortality.  PDF of description is included. 

# In this script I am seeking to replicate Tom's work.  Next steps include the
# possibility of a spline, or a bayesian application.  I will also run this
# estimator on as many curves as possible. 

# Note that if any of the estimators are outside of the bounds of data,
# the function will break.  I do not have a good solution for this as of yet. 
#------------------------------------------------------------------------------#

##### Load data #####
library(tidyverse)
library(magrittr)
library(numDeriv)
setwd("~/S-timators")
rm(list=ls())
# gc()
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

# ggplot() + 
#   geom_line(data = m.sim, aes(x = age, y = v), color = "blue") + 
#   geom_line(data = f.sim, aes(x = age, y = v), color = "red")  

# It can be seen that the female curve is right-shifted and scaled
# down relative to the male curve


#------------------------------------------------------------------------------#

procedure <- function(m.sim, f.sim, ALPHA, BETA, n = 10) {
  ##### Step 0: Normalizing densities #####
  # I made this mistake before by not normalizing.  
  m.sim$v <- m.sim$v/sum(m.sim$v)
  f.sim$v <- f.sim$v/sum(f.sim$v)
  
  ##### Implementation #####
  # Notation may be slightly different from the pdf
  # all functions have prefix f_{quantity}
  # CONSTANTS are in CAPITALS
  # quantities should be consistent like: {variable}.{superscript}_{sex} sigma.sq.hat_m
  
  
  #### Convert to continuous time functions #
  v_w <- splinefun(m.sim$age, m.sim$v)
  v_m <- splinefun(f.sim$age, f.sim$v)
  # v_w <- approxfun(m.sim$age, m.sim$v)
  # v_m <- approxfun(f.sim$age, f.sim$v)
  
  # ALPHA <- min(min(m.sim$age), min(f.sim$age))
  # BETA <- max(max(m.sim$age), max(f.sim$age))

  m.sim <- m.sim %>% filter(age >= ALPHA, age <= BETA)
  f.sim <- f.sim %>% filter(age >= ALPHA, age <= BETA)
  mu_m <- weighted.mean(m.sim$age, m.sim$v)
  mu_f <- weighted.mean(f.sim$age, f.sim$v)
  
  DELTA <- round(mu_f - mu_m)
  
  #### Iteration ####
  # This function is a wrapper that returns all the constants needed for the
  # estimators.  The default value of DELTA is mu_f - mu_m. 
  
  
  
  stimator <- function(DELTA, ALPHA, BETA, v_w, v_m){
    
    if (is.nan(DELTA)) {message("DELTA NAN"); return(0)}
    if (ALPHA + DELTA > BETA) {message("ALPHA + DELTA > BETA"); return(0)}
    if (BETA - DELTA < ALPHA) {message("BETA - DELTA < ALPHA"); return(0)}
    
    #### mu ####
    
    mu.hat_w <- integrate(function(a) {a * v_w(a)}, ALPHA+DELTA, BETA)$value
    mu.hat_m <- integrate(function(a) {a * v_m(a)}, ALPHA, BETA-DELTA)$value
    
    # Quick sanity check:
    # ggplot() + 
    #   geom_line(data = m.sim, aes(x = age, y = v), color = "blue") + 
    #   geom_line(data = f.sim, aes(x = age, y = v), color = "red")  +
    #   geom_vline(xintercept = ALPHA + DELTA) + 
    #   geom_vline(xintercept = BETA - DELTA) + 
    #   geom_vline(xintercept = mu_m, color = "blue") + 
    #   geom_vline(xintercept = mu_f, color = "red") 
    
    #### sigma.hat.sq ####
    
    tryCatch({
      sigma.hat.sq_w <- integrate(function(a){(a - mu.hat_w)^2 * v_w(a)},ALPHA+DELTA, BETA)$value
      sigma.hat.sq_m <- integrate(function(a){(a - mu.hat_m)^2 * v_m(a)},ALPHA, BETA-DELTA)$value
    }, error = function(e) {return(0)})
    
    
    #### Q and Q.tilde ####
    
    Q <- integrate(v_w,ALPHA + DELTA, BETA)$value
    Q.tilde <- integrate(v_m, ALPHA, BETA - DELTA)$value
    
    #### R and R.tilde ####
    
    R <- mu.hat_m + DELTA
    R.tilde <- mu.hat_w - DELTA
    
    
    #### N and N.tilde ####
    
    N <- ((BETA - R)^2 - sigma.hat.sq_m/Q.tilde)*v_w(BETA) - 
      ((ALPHA - mu.hat_m)^2 - sigma.hat.sq_m/Q.tilde)*v_w(ALPHA + DELTA) - 
      2*mu.hat_w+2*R*Q
    
    N.tilde <- ((BETA - mu.hat_w)^2 - sigma.hat.sq_w/Q)*v_m(BETA - DELTA) - 
      ((ALPHA - R.tilde)^2 - sigma.hat.sq_w/Q)*v_m(ALPHA) -
      2*mu.hat_m + 2*R.tilde*Q.tilde
    
    #### D and D.tilde ####
    # Derivative is v_m(x, deriv = 1)
    
    D <- ((BETA - R)^2 - sigma.hat.sq_m/Q.tilde)*v_m(DELTA, deriv = 1) - 
      ((ALPHA - mu.hat_m)^2-sigma.hat.sq_m/Q.tilde)*v_m(ALPHA + DELTA, deriv = 1) -
      2*(BETA - R)*v_w(BETA) + 2*(ALPHA - mu.hat_m)*v_w(ALPHA + DELTA) + 2*Q
    
    
    D.tilde <- ((BETA - mu.hat_w)^2 - sigma.hat.sq_w/Q)*v_m(BETA - DELTA, deriv = 1)-
      ((ALPHA - R.tilde)^2 - sigma.hat.sq_w/Q)*v_m(ALPHA, deriv = 1) - 
      2*(BETA - mu.hat_w)*v_m(BETA - DELTA) + 2*(ALPHA - R.tilde) * v_m(ALPHA)+2*Q.tilde
    
    
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
  
  ret <- stimator(DELTA, ALPHA, BETA, v_m, v_w)
  
  #### The Formulae  ####
  #### Estimator 1 ####
  f_e1 <- function(delta, 
                   N, 
                   D, 
                   sigma.hat.sq_w , 
                   sigma.hat.sq_m , 
                   mu.hat_w , 
                   Q , 
                   Q.tilde , 
                   R) {
    
    a1 <- delta - N/D
    a2 <- N^2 - 2*D*(sigma.hat.sq_w + 2*mu.hat_w^2 - Q*mu.hat_w^2 - 2*R*mu.hat_w + 
                       Q*R^2 -  Q/Q.tilde*sigma.hat.sq_m)
    return(a1 - (sqrt(a2)/D))
  }
  
  #### Estimator 2 ####
  f_e2 <- function(delta, 
                   N.tilde ,
                   D.tilde ,
                   sigma.hat.sq_m ,
                   mu.hat_m ,
                   Q.tilde ,
                   R.tilde ,
                   Q ,
                   sigma.hat.sq_w) {
    
    a1 <- delta + N.tilde/D.tilde
    a2 <-  N.tilde^2 - 2*D.tilde * (sigma.hat.sq_m + 2*mu.hat_m^2 - 
                                      Q.tilde*mu.hat_m^2 - 2*R.tilde*mu.hat_m +
                                      Q.tilde * R.tilde^2 - Q.tilde/Q*sigma.hat.sq_w)
    
    return(a1 + sqrt(a2)/D.tilde)
    
  }
  
  #### Estimator 3 ####
  
  f_e3 <- function(delta,
                   sigma.hat.sq_w ,
                   mu.hat_w ,
                   Q ,
                   R , 
                   Q.tilde ,
                   sigma.hat.sq_m ,
                   N) {
    
    a1 <- sigma.hat.sq_w + 2*mu.hat_w^2 - Q*mu.hat_w^2 - 2*R*mu.hat_w +
      Q*R^2 - Q/Q.tilde*sigma.hat.sq_m
    # browser()
    return(delta - (a1/N))
  }
  
  #### Estimator 4 ####
  
  f_e4 <- function(delta ,
                   sigma.hat.sq_m ,
                   mu.hat_m,
                   Q.tilde ,
                   R.tilde ,
                   Q ,
                   sigma.hat.sq_w ,
                   N.tilde ) {
    a1 <- sigma.hat.sq_m + 2*mu.hat_m^2 - Q.tilde*mu.hat_m^2 - 2*R.tilde*mu.hat_m +
      Q.tilde * R.tilde^2 - Q.tilde/Q*sigma.hat.sq_w
    # browser()
    return(delta + (a1/N.tilde))
  }
  
  
  #### Iteration Procedure ####
  # Program has already run once, providing first estimates of e1, e2, e3, e4
  # Create nx4 matrix of values
  # Call iteratively replacing DELTA with e1
  # After completing n iterations, run to clear, and repeat with e2, e3, e4
  
  iteration_df <- data.frame(matrix(NA, n, 5))
  colnames(iteration_df) <- c("index", "e1", "e2", "e3", "e4")
  iteration_df[1, 2:5] <- DELTA
  iteration_df$index <- 1:n
  
  # e1
  for (i in 2:n) {
    ret <- stimator(DELTA = iteration_df[i - 1, "e1"], ALPHA, BETA, v_m, v_w)
    if (length(ret) == 1) {break}
    e1 <- with(ret, f_e1(delta = iteration_df[i - 1, "e1"],
                         N, 
                         D, 
                         sigma.hat.sq_w , 
                         sigma.hat.sq_m , 
                         mu.hat_w , 
                         Q , 
                         Q.tilde , 
                         R))
    iteration_df[i, "e1"] <- e1
  }
  
  # e2 
  for (i in 2:n) {
    ret <- stimator(DELTA = iteration_df[i - 1, "e2"], ALPHA, BETA, v_m, v_w)
    if(length(ret) == 1) {break}
    e2 <- with(ret, f_e2(delta = iteration_df[i - 1, "e2"], 
                         N.tilde = N.tilde,
                         D.tilde = D.tilde,
                         sigma.hat.sq_m = sigma.hat.sq_m,
                         mu.hat_m = mu.hat_m,
                         Q.tilde = Q.tilde,
                         R.tilde = R.tilde,
                         Q = Q,
                         sigma.hat.sq_w = sigma.hat.sq_w))
    iteration_df[i, "e2"] <- e2
  }
  
  # e3
  for (i in 2:n) {
    ret <- stimator(DELTA = iteration_df[i - 1, "e3"], ALPHA, BETA, v_m, v_w)
    if (length(ret) == 1 ) {break}
    e3 <- with(ret, f_e3(delta = DELTA,
                         sigma.hat.sq_w = sigma.hat.sq_w,
                         mu.hat_w = mu.hat_w,
                         Q = Q,
                         R = R, 
                         Q.tilde = Q.tilde,
                         sigma.hat.sq_m = sigma.hat.sq_m,
                         N=N))
    iteration_df[i, "e3"] <- e3
  }
  
  # e4
  for (i in 2:n) {
    ret <- stimator(DELTA = iteration_df[i-1, "e4"], ALPHA, BETA, v_m, v_w)
    if (length(ret) == 1) {break }
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
  
  return(iteration_df)
}

(procedure(m.sim, f.sim, ALPHA = 70, BETA = 90, n = 10))
(procedure(m.sim, f.sim, 70, 80, 10))
(procedure(m.sim, f.sim, 70, 90, 10))
(procedure(m.sim, f.sim, 80, 90, 10))
(procedure(m.sim, f.sim, 80, 100, 10))
(procedure(m.sim, f.sim, 50, 100, 100))
(procedure(m.sim, f.sim, 30, 95, 10))

# Create a contour plot with 5-year resolution

contour_df <- data.frame(ALPHA = NA, 
                         BETA  = NA,
                         e1 = NA,
                         e2 = NA, 
                         e3 = NA, 
                         e4 = NA)

for (A in seq(0, 120, 5)){
  for (B in seq(0, 120, 5)){
    if (B <= A) {next}
    # print(c(A, B))
    t <- suppressMessages(suppressMessages(procedure(m.sim, f.sim, A, B, n = 100)))
    t <- t %>% slice(100) %>% select(e1, e2, e3, e4)
    t <- t %>% mutate(ALPHA = A, BETA = B)
    contour_df <- contour_df %>% bind_rows(t)
  }
}


f <- scale_fill_gradient2(low = "red", mid = "white",
                     high = "blue", midpoint = 7, space = "Lab",
                     limits = c(-10, 10), 
                     na.value = "grey50", guide = "colourbar", aesthetics = "fill")

gridExtra::grid.arrange(
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e1)) +
    f +
    ggtitle("e1"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e2)) + 
    f + 
    ggtitle("e2"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e3)) + 
    f +
    ggtitle("e3"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e4)) +
    f + 
    ggtitle("e4")
) + ggtitle("100 Iteration Contour Plot (Values >10, <-10 removed)")




SWE <- read_table("SWE_mx.txt", skip = 1)

# Some data reformating
SWE <- SWE %>% mutate(Age = if_else(Age== "110+", "110", Age),
                      Age = as.numeric(Age),
                      Male = if_else(Male == ".", as.character(NA), Male),
                      Male = as.numeric(Male),
                      Female = if_else(Female == ".", as.character(NA), Female),
                      Female = as.numeric(Female),
                      Total = if_else(Total == ".", as.character(NA), Total),
                      Total = as.numeric(Total))

# Drop non-complete cohorts
SWE <- SWE %>% group_by(Year) %>%
  filter(!any(is.na(Male))) %>%
  filter(!any(is.na(Female))) %>%
  filter(!any(is.na(Total)))  %>%
  ungroup()

# Apparently this only gives the 1892 cohort, which is fine for now

m.swe <- SWE %>% select(Age, Male) %>%
  rename(age = Age, v = Male) 
m.swe  

f.swe <- SWE %>% select(Age, Female) %>%
  rename(age = Age, v = Female)
f.swe

procedure(m.swe, f.swe, 60, 100, 10)  # e1 converges around 11 but e4 to 33
procedure(m.swe, f.swe, 80, 100, 100) #e1 converges to 13 but e2 converes to -8?
procedure(m.swe, f.swe, 80, 90, 100) # e1 converges to 5 but e3 converges to 3.84


# Create a contour plot with 5-year resolution

contour_df <- data.frame(ALPHA = NA, 
                         BETA  = NA,
                         e1 = NA,
                         e2 = NA, 
                         e3 = NA, 
                         e4 = NA)

for (A in seq(0, 120, 5)){
  for (B in seq(0, 120, 5)){
    if (B <= A) {next}
    # print(c(A, B))
    t <- suppressMessages(suppressMessages(procedure(m.swe, f.swe, A, B, n = 100)))
    t <- t %>% slice(100) %>% select(e1, e2, e3, e4)
    t <- t %>% mutate(ALPHA = A, BETA = B)
    contour_df <- contour_df %>% bind_rows(t)
  }
}


f <- scale_fill_gradient2(low = "red", mid = "white",
                          high = "blue", midpoint = 0, space = "Lab",
                          limits = c(-10, 10), 
                          na.value = "grey50", guide = "colourbar", aesthetics = "fill")

gridExtra::grid.arrange(
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e1)) +
    f +
    ggtitle("e1"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e2)) + 
    f + 
    ggtitle("e2"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e3)) + 
    f +
    ggtitle("e3"),
  contour_df %>% 
    ggplot() + 
    geom_tile(mapping = aes(x = ALPHA, y = BETA, fill = e4)) +
    f + 
    ggtitle("e4")
) + ggtitle("100 Iteration Contour Plot (Values >10, <-10 removed)")







  