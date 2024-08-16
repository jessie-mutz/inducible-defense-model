## Mutz, J. & K.C. Abbott. Life history modulates effects of inducible defenses 
##  on consumer-resource dynamics.
## Supplemental code: Simulations
## JM, Last updated 2024-06-05

# Set working directory.
#setwd()

# Source model functions.
source("./scripts/01_model-functions.R")

## Set up dataframe for combinations of plant maturation rate (g), reproductive 
## plant mortality (m), and per-capita plant fecundity (b) to represent plant
## populations across a range of life histories.================================

# Generate 500 samples of plant maturation rate (g) and reproductive plant mortality
#  rate (m) for simulations.
set.seed(1408)
samples <- data.frame(cbind(
  g = runif(500, min = 0.05, max = 0.95),
  m = runif(500, min = 0.05, max = 0.95)
))

# Write function to calculate threshold value of per-capita plant fecundity (b)
#  necessary for population persistence in the absence of inducible defenses,
#  for given values of g and m
b_func <- function(g, m, w = 1, u = 0.15, x = 1) {
  b <- (1/(g*x)) * (w - u + g) * (w - u + m)
  return(b)
}

# Use 'b_func' to assign b values for each sampled set of m and g values.
samples$b <- b_func(g = samples$g, m = samples$m)

## Run simulations to investigate the effects of inducible defense costs (theta),
## inducible defense responsiveness (tau), and inducible defense effectiveness
## (epsilon) across plant life histories.=======================================

# Function for running simulations: for each element of 'samples': 
#  reassign b, m, and g; use 'funcy' to find equilibrium as a function of
#  variables 'funcy_a', 'funcy_b', and 'funcy_c'; and return qualitative and 
#  quantitative outputs as separate objects.

sim_func <- function(funcy_a, funcy_b, funcy_c) {
  
  samples <- get('samples')
  params <- get('params')
  param.vals <- get('param.vals')
  qual_out <- quan_out <- list(NA)
  times <- seq(0, 5000, by = 1)
  
  for(i in 1:length(samples[,1])) {
    
    # Reassign b, m, and g according to 'samples' df
    params[c('b', 'm', 'g')] <- c(samples$b[i], samples$m[i], samples$g[i])
    
    # Use 'funcy' to find equilibrium as a function of theta, tau, and epsilon values
    #  specified in the 'param.vals' list
    out1 <- funcy(a = funcy_a, b = funcy_b, c = funcy_c, temp.params = params)
    
    
    # Melt qualitative output array, add replicate number, and add parameter values
    out.qual.melt <- melt(out1[[2]])
    out.qual2 <- reshape2::dcast(data = out.qual.melt, 
                                 formula = out.qual.melt[,1] + out.qual.melt[,2] + 
                                   out.qual.melt[,3] ~ qual, 
                                 margin = "value")
    colnames(out.qual2)[1:3] <- c(funcy_a[1], funcy_b[1], funcy_c[1])
    out.qual3 <- data.frame(cbind(
      rep = i,
      out1[[3]],
      out.qual2
    ))
    
    # Melt quantitative output array, add replicate number, and add parameter values
    out.quan.melt <- melt(out1[[1]])
    out.quan2 <- reshape2::dcast(data = out.quan.melt, 
                                 formula = out.quan.melt[,1] + out.quan.melt[,2] + 
                                   out.quan.melt[,3] + timestep ~ variable, 
                                 margin = "value")
    colnames(out.quan2)[1:3] <- c(funcy_a[1], funcy_b[1], funcy_c[1])
    out.quan3 <- data.frame(cbind(
      rep = i,
      out1[[3]],
      out.quan2
    ))
    
    # Add columns for prop. of mature plant biomass (P.prop), mean induction of
    #  seedlings (i_s_mean) and mature plants (i_p_mean), herbivore density per
    #  unit plant biomass (H_SP), and total no. of seeds produced (seeds).
    out.quan3[c('P.prop', 'i_s_mean', 'i_p_mean', 'H_SP', 'seeds')] <- c( 
      with(out.quan3, P / (S + P)),
      with(out.quan3, I_s / S),
      with(out.quan3, I_p / P),
      with(out.quan3, H / (S + P)),
      with(out.quan3, b*(P - I_p*theta)))
    out.quan3$seeds <- sapply( out.quan3$seeds, function(x) { max(0, x) } )
    
    # Save each rep as a separate list element
    qual_out[[i]] <- out.qual3
    quan_out[[i]] <- out.quan3
    
    cat(i, "\n") 
  }
  
  return(list(qual_out, quan_out)) 
  
} # end of function

# Define initial parameter list.
params <- list(
  b = 6,
  theta = 1.5,
  g = 0.15,
  x = 1,
  tau_s = 0.5,
  tau_p = 0.5,
  c = 0.8,
  a = 0.1,
  ep_s = 30,
  ep_p = 30,
  w = 1,
  r = 2,
  mu = 0.15,
  m = 0.05
)

## Run simulations.=============================================================

###### NOT RUN (provided in "output/01_tau-theta-ep27-qual.csv"):
## Simulate across values of inducible defense responsiveness (tau) and inducible 
##  defense costs (theta) for 500 plant life history scenarios (i.e., combination 
##  of seedling maturation rate [g] and reproductive plant mortality rate [m]).
#param.vals <- list(
#  'theta' = seq(2, 4, by = 0.5),
#  'tau_s' = seq(0.2, 1, by = 0.2),
#  'tau_p' = seq(0.2, 1, by = 0.2),
#  'ep_s' = 27,
#  'ep_p' = 27
#)
## Run simulations across parameter values specified above.
#out <- sim_func( 'theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
#qual_out <- out[[1]] #assign qualitative simulation output to its own object
## Convert qualitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#qual_df <- do.call(rbind.data.frame, qual_out)
#write.csv(qual_df, file = "output/01_tau-theta-ep27-qual.csv", row.names = FALSE) #save dataframe as .csv
#rm(qual_out, qual_df) #remove objects from environment

###### NOT RUN (provided in "output/02_tau0.8-theta-ep-qual.csv"):
## Simulate across values of inducible defense costs (theta) and inducible defense 
##  efficacy (epsilon) for 500 plant life history scenarios (i.e., combination 
##  of seedling maturation rate [g] and reproductive plant mortality rate [m]).
#param.vals <- list(
#  'theta' = c(2, 4),
#  'tau_s' = 0.8,
#  'tau_p' = 0.8,
#  'ep_s' = seq(5, 45),
#  'ep_p' = seq(5, 45)
#)
## Run simulations across parameter values specified above.
#out <- sim_func( 'theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
#qual_out <- out[[1]] #assign qualitative simulation output to its own object
## Convert qualitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#qual_df <- do.call(rbind.data.frame, qual_out)
#write.csv(qual_df, file = "output/02_tau0.8-theta-ep-qual.csv", row.names = FALSE) #save dataframe as .csv
#rm(qual_out, qual_df) #remove objects from environment

###### NOT RUN (provided in "output/03_tau0.8-theta4-ep-g0.955-qual.csv"):
## Simulate across values of inducible defense efficacy (epsilon) at increased 
##  seedling maturation rate (g = 0.955) for 500 plant life history scenarios.
#samples$g <- 0.955 #reassign values of g to increased value (0.955)
#param.vals <- list(
#  'theta' = 4,
#  'tau_s' = 0.8,
#  'tau_p' = 0.8,
#  'ep_s' = seq(5, 45),
#  'ep_p' = seq(5, 45)
#)
## Run simulations across parameter values specified above.
#out <- sim_func( 'theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
#qual_out <- out[[1]] #assign qualitative simulation output to its own object
## Convert qualitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#qual_df <- do.call(rbind.data.frame, qual_out)
#write.csv(qual_df, file = "output/03_tau0.8-theta4-ep-g0.955-qual.csv", row.names = FALSE) #save dataframe as .csv
#rm(qual_out, qual_df) #remove objects from environment

###### NOT RUN (provided in "output/04_tau-theta-ep27-r8-qual.csv"):
## Simulate across values of inducible defense responsiveness (tau) and inducible 
##  defense costs (theta) at high herbivore intrinsic growth rate (r = 8) for 
##  500 plant life history scenarios (i.e., combination of seedling maturation 
##  rate [g] and reproductive plant mortality rate [m]).
#params['r'] <- 8
#param.vals <- list(
#  'theta' = seq(2, 4, by = 0.5),
#  'tau_s' = seq(0.2, 1, by = 0.2),
#  'tau_p' = seq(0.2, 1, by = 0.2),
#  'ep_s' = 27,
#  'ep_p' = 27
#)
## Run simulations across parameter values specified above.
#out <- sim_func( 'theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
#qual_out <- out[[1]] #assign qualitative simulation output to its own object
## Convert qualitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#qual_df <- do.call(rbind.data.frame, qual_out)
#write.csv(qual_df, file = "output/04_tau-theta-ep27-r8-qual.csv", row.names = FALSE) #save dataframe as .csv
#rm(qual_out, qual_df) #remove objects from environment

## For quantitative simulations, use nine representative combinations of seedling
##  maturation rate (g) and reproductive plant mortality (m):
samples <- data.frame(expand.grid(
  g = c(0.15, 0.5, 0.8),
  m = c(0.05, 0.5, 0.8)
))

# Use 'b_func' to assign b values for each sampled set of m and g values.
samples$b <- b_func(g = samples$g, m = samples$m)

###### NOT RUN (provided in "output/05_tau-theta-ep-quan.csv"):
## Simulate across values of inducible defense responsiveness (tau), inducible 
##  defense costs (theta), and inducible defense efficacy (epsilon) for 9 
##  representative plant life history scenarios (i.e., combination of seedling 
##  maturation rate [g] and reproductive plant mortality rate [m]).
#param.vals <- list(
#  'theta' = c(2, 4),
#  'tau_s' = seq(0, 1, by = 0.05),
#  'tau_p' = seq(0, 1, by = 0.05),
#  'ep_s' = c(15, 30, 45),
#  'ep_p' = c(15, 30, 45)
#)
## Run simulations across parameter values specified above.
#out <- sim_func( 'theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
#quan_out <- out[[2]] #assign quantitative simulation output to its own object
## Convert quantitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#quan_df <- do.call(rbind.data.frame, quan_out)
#write.csv(quan_df, file = "output/05_tau-theta-ep-quan.csv", row.names = FALSE) #save dataframe as .csv
#rm(quan_out, quan_df) #remove objects from environment

###### NOT RUN (provided in "output/06_taus-taup-theta1.5-ep-quan.csv"):
## Simulate across values of inducible defense responsiveness of seedlings (tau_s) 
##  and reproductive plants (tau_p) for 4 representative plant life history scenarios 
##  (i.e., combination of seedling maturation rate [g] and reproductive plant 
##  mortality rate [m]).
#param.vals <- list(
#  'ep_s' = 30,
#  'ep_p' = 30,
#  'tau_s' = seq(0, 1, by = 0.05),
#  'tau_p' = seq(0, 1, by = 0.05)
#)
## Run simulations across parameter values specified above.
#out <- sim_func( c('ep_s', 'ep_p'), 'tau_s', 'tau_p')
#quan_out <- out[[2]] #assign quantitative simulation output to its own object
## Convert quantitative simulation output from large list into dataframe with 
## column 'rep' designating list element.
#quan_df <- do.call(rbind.data.frame, quan_out)
#quan_df <- subset(quan_df, timestep == 100)
#write.csv(quan_df, file = "output/06_taus-taup-theta1.5-ep30-quan.csv", row.names = FALSE) #save dataframe as .csv
#rm(quan_out, quan_df) #remove objects from environment
