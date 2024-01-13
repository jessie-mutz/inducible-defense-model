## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Qualitative simulations
## JM, Last updated 2024-01-13

# Source model functions for model simulations.
source("./inducible-defense-model-v0/01_model-functions.R")

# Because the same set of qualitative simulations is desired for multiple plant 
#  demography scenarios, specify a function for simulating:
#   (A) Qualitative behavior ~ tau x epsilon,
#   (B) Qualitative behavior ~ tau x b,
#   (C) Qualitative behavior ~ tau x g, 
#   (D) Qualitative behavior ~ tau x m,
#  each at baseline and increased theta, increased c, and increased r.
# Output from each of the four simulation sets is saved in a list, and the four 
#  lists are saved as separate RData files (filename provided as an argument).
# Arguments:
#  params = vector of named parameter values;
#  yini = initial values for state variables;
#  filename = character string for name of list to save;
#  n.times = number of timesteps;
#  tau.range = vector of tau values to use;
#  ep.range = vector of epsilon values to use;
#  b.range = vector of b values to use;
#  g.range = vector of g values to use;
#  m.range = vector of m values to use.

qual_set <- function(params, yini, filename, n.times = 5000,
                     tau.range = seq(0, 1, by = 0.02), 
                     ep.range = seq(0, 30, by = 1), 
                     b.range = seq(0, 15, by = 0.5), 
                     g.range = seq(0, 1, by = 0.02), 
                     m.range = seq(0.02, 1, by = 0.02)) {
  
  tau_ep_out <- tau_b_out <- tau_g_out <- tau_m_out <- list(NA)
  params.temp <- params
  
  # Simulations at baseline and increased theta.
  thetas <- c(1.5, 4)
  dim.temp <- 1
  for( i in 1:2 ) {
    params.temp[c('theta')] <- thetas[i]
    tau_ep_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                       x.axis = 'tau', x.range = tau.range, 
                                       y.axis = 'ep', y.range = ep.range)
    tau_b_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                      x.axis = 'tau', x.range = tau.range, 
                                      y.axis = 'b', y.range = b.range, exactY = TRUE)
    tau_g_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                      x.axis = 'tau', x.range = tau.range, 
                                      y.axis = 'g', y.range = g.range, exactY = TRUE)
    tau_m_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                      x.axis = 'tau', x.range = tau.range, 
                                      y.axis = 'm', y.range = m.range, exactY = TRUE)
    dim.temp <- dim.temp + 1
  }
  
  # Simulations at increased c.
  params.temp <- params #reset params
  params.temp[c('c')] <- 5
  tau_ep_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                     x.axis = 'tau', x.range = tau.range, 
                                     y.axis = 'ep', y.range = ep.range)
  tau_b_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'b', y.range = b.range, exactY = TRUE)
  tau_g_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'g', y.range = g.range, exactY = TRUE)
  tau_m_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'm', y.range = m.range, exactY = TRUE)
  dim.temp <- dim.temp + 1
  
  # Simulations at increased r.
  params.temp <- params #reset params
  params.temp[c('r')] <- 8
  tau_ep_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                     x.axis = 'tau', x.range = tau.range, 
                                     y.axis = 'ep', y.range = ep.range)
  tau_b_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'b', y.range = b.range, exactY = TRUE)
  tau_g_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'g', y.range = g.range, exactY = TRUE)
  tau_m_out[[dim.temp]] <- qual_out(n.timesteps = n.times, param.vec = params.temp, 
                                    x.axis = 'tau', x.range = tau.range, 
                                    y.axis = 'm', y.range = m.range, exactY = TRUE)
  # Name list elements accordingly.
  names(tau_ep_out) <- c("baseline", "theta = 4", "c = 5", "r = 8")
  names(tau_b_out) <- c("baseline", "theta = 4", "c = 5", "r = 8")
  names(tau_g_out) <- c("baseline", "theta = 4", "c = 5", "r = 8")
  names(tau_m_out) <- c("baseline", "theta = 4", "c = 5", "r = 8")
  
  # Save all lists as RData.
  save(tau_ep_out, tau_b_out, tau_g_out, tau_m_out, file = filename)
  
} #end function

# Initial parameter list.
params <- c(
    b = 6,
    theta = 1.5,
    g = 0.15,
    x = 1,
    tau_s = 0.5,
    tau_p = 0.5,
    c = 0.8,
    a = 0.1,
    ep_s = 25,
    ep_p = 25,
    w = 1,
    r = 2,
    mu = 0.15,
    m = 0.05
  )

# Set initial values for state variables.
yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001)

## Scenario 1: Low mature plant mortality/low seedling maturation rate.=========
params[c('b','g','m')] <- c(6, 0.15, 0.05)

# Run 'qual_set' function for this set of parameter values.
## NOT RUN --> output provided in '07_qualitative-output-lowm-lowg.RData' ##
#qual_set(params = params, yini = yini, 
#         filename = paste0(Sys.Date(), "_qualitative-output-lowm-lowg.RData"))

## Scenario 2: High mature plant mortality/low seedling maturation rate.========
params[c('b','g','m')] <- c(11, 0.15, 0.8)

# Run 'qual_set' function for this set of parameter values.
## NOT RUN --> output provided in '08_qualitative-output-highm-lowg.RData' ##
#qual_set(params = params, yini = yini, 
#         filename = paste0(Sys.Date(), "_qualitative-output-highm-lowg.RData"))

## Scenario 3: Low mature plant mortality/moderate seedling maturation rate.====
params[c('b','g','m')] <- c(2.43, 0.5, 0.05)

# Run 'qual_set' function for this set of parameter values.
## NOT RUN --> output provided in '09_qualitative-output-lowm-modg.RData' ##
#qual_set(params = params, yini = yini, 
#         filename = paste0(Sys.Date(), "_qualitative-output-lowm-modg.RData"))

## Scenario 4: High mature plant mortality/moderate seedling maturation rate.===
params[c('b','g','m')] <- c(4.445, 0.5, 0.8)

# Run 'qual_set' function for this set of parameter values.
## NOT RUN --> output provided in '10_qualitative-output-highm-modg.RData' ##
#qual_set(params = params, yini = yini, 
#         filename = paste0(Sys.Date(), "_qualitative-output-highm-modg.RData"))
