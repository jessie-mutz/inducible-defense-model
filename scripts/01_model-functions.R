## Mutz, J. & K.C. Abbott. Life history modulates effects of inducible defenses 
##  on consumer-resource dynamics.
## Supplemental code: Model functions
## JM, Last updated 2024-06-05

# Load required packages.
require(deSolve)
require(parallel)
require(reshape2)

# Specify ODE.
full <- function (Time, State, Params) {
  with( as.list( c(State, Params) ), {
    dam_s <- w * H / ( (S + P) * ( 1 + ep_s * (I_s/S) ) + H )
    dam_p <- w * H / ( (S + P) * ( 1 + ep_p * (I_p/P) ) + H )
    dam.per.H <- (dam_s * S + dam_p * P) / H
    Seedlings	<- max(0, b * (P - theta * I_p) * (1 - P)) - g * S - dam_s * S
    Mature.plants <- g * x * S - dam_p * P - m * P
    Induction_S <- tau_s * max(0, (S - I_s)) * ( dam_s / (c + dam_s) ) - 
      a * I_s - g * I_s - dam_s * I_s
    Induction_P <- tau_p * max(0, (P - I_p)) * ( dam_p / (c + dam_p) ) - 
      a * I_p - m * I_p + g * x * ( ep_s / ep_p ) * I_s - dam_p * I_p
    Herbivores <- r * H * ( (dam.per.H^2 - mu^2) / (dam.per.H^2 + mu^2) )
    
    return(list(c(Seedlings, Mature.plants, Induction_S, Induction_P, Herbivores)))
  } )
}

## Functions for determining qualitative behaviors and quantifying densities at
## equilibrium. ===============================================================
# For each combination of parameter values, run the ODE across many timesteps and
#  (1) use the last 10% of timesteps to characterize qualitative outcomes for 
#      each state variable: positive (TRUE/FALSE), stable equilibrium 
#      (TRUE/FALSE), cyclic dynamics (according to 2 different diagnostics),
#  (2) quantify equilibrium densities for each state variable.
# Requires defined objects 'times' (vector of timesteps), 'params' (vector of
#  named parameter values), and 'param.vals' (list including vectors named a, b, 
#  and c, each containing range of parameter values).

funcy <- function(a, b, c, temp.params = params) {
  
  # Set up arrays for qualitative and quantitative outcomes.
  times <- get('times')
  param.vals <- get('param.vals')
  
  length.out <- round(0.02*length(times))
  
  out.array <- array( NA, dim = c( length( param.vals[[a[1]]] ),
                                   length( param.vals[[b[1]]] ), 
                                   length( param.vals[[c[1]]] ), 
                                   length.out, 5) )
  out.qual <- array( NA, dim = c(dim(out.array)[1:3], 31))
  
  # Name columns for keeping track of qualitative behavior of each state variable 
  # at equilibrium.
  qual.names <- c(sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".pos") } ),
                  sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".neg") } ),
                  sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".pt.eq") } ),
                  sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".cyc1") } ),
                  sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".cyc2") } ),
                  sapply( c('S', 'P', 'I_s', 'I_p', 'H'), 
                          FUN = function(x) { paste0(x, ".NA") } ),
                  "qual")
  
  all.names <- list( param.vals[[a[1]]], param.vals[[b[1]]], param.vals[[c[1]]], 
                     1:length.out, c('S', 'P', 'I_s', 'I_p', 'H'),
                     qual.names)
  names(all.names) <- list( a[1], b[1], c[1], 'timestep', 'variable', 'qual' )
  dimnames(out.array) <- all.names[1:5]
  dimnames(out.qual) <- all.names[c(1:3,6)]
  
  
  # Solve ODE and characterize behavior and pattern of densities at equilibrium.
  for ( x in 1:length(param.vals[[a[1]]]) ) {
    for ( y in 1:length(param.vals[[b[1]]]) ) {
      for ( z in 1:length(param.vals[[c[1]]]) ) {
        
        # Set parameter values using values of a, b, and c.
#        temp.params <- get('params')
        temp.params[a] <- param.vals[[a[1]]][x]
        temp.params[b] <- param.vals[[b[1]]][y]
        temp.params[c] <- param.vals[[c[1]]][z]
        
        skip_to_next <- FALSE
        
        tryCatch({
          
          # Set initial values for state variables and use ODE solver to get
          #  densities through time.
          yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001)
          out <- deSolve::ode(y = yini, times = times, func = full, parms = temp.params)
          to.length <- min(dim(out)[1], length.out)
          
          # Save equilibrium densities to 'out.array' (or skip if error encountered).
          out.array[x,y,z,,] <- out[(dim(out)[1]-(to.length-1)):(dim(out)[1]), 
                                    c('S', 'P', 'I_s', 'I_p', 'H')]
          
          # Focus on the last 10% of timesteps to characterize qualitative behavior.
          to.length2 <- min(dim(out)[1], round(0.1*length(times)))
          last.bit <- data.frame(out[(dim(out)[1]-(to.length2-1)):(dim(out)[1]), 
                                     c('S', 'P', 'I_s', 'I_p', 'H')])
          
          # Characterize qualitative outcomes for each state variable: 
          #  positive (TRUE/FALSE), stable equilibrium (TRUE/FALSE), cyclic 
          #  dynamics (according to 2 different diagnostics).
          
          var_pos <- sapply(last.bit[,labels(yini)], function(x) { all( x > 1e-08 ) })
          var_neg <- sapply(last.bit[,labels(yini)], function(x) { any( x < 0 ) })
          var_pt.eq <- sapply(last.bit[,labels(yini)], function(x) { isTRUE( max(x) - min(x) < 1e-03 ) })
          var_cyc1 <- sapply(last.bit[,labels(yini)], function(x) { 
            isTRUE( abs(which(x == min(x)) - which(x == max(x))) < length(last.bit[,1]) - 1) })
          var_cyc2 <- sapply(last.bit[,labels(yini)], function(x) { sd(x) / mean(x) })
          var_NA <- sapply(last.bit[,labels(yini)], function(x) { any( is.na(x) == TRUE ) })
          
          # Summarize qualitative outcomes: positive point equilibrium (all state 
          #  variables), extinction (all state variables), cycles (all state 
          #  variables), or NA (none of the above, including differences in 
          #  qualitative outcomes among state variables).
          qual <- NA
          if (sum(var_NA) == 5) { qual <- "unk" }
          if (sum(var_pos) == 5 & sum(var_pt.eq) == 5) { qual <- "pos.pt.eq" }
          if (sum(var_pos) == 0) { qual <- "extinct" }
          if (sum(var_cyc1) == 5 & sum(var_pos) == 5) { qual <- "cycles" }
          
          out.qual[x,y,z,] <- c(var_pos, var_neg, var_pt.eq, var_cyc1, var_cyc2, var_NA, qual)
        }, error = function(e) { skip_to_next <<- TRUE })
        
        if(skip_to_next) { next }  
        
        # Print parameter values used in current run of loop to track function progress.
        cat(paste(a, " = ", param.vals[[a[1]]][x]), 
            paste(b, " = ", param.vals[[b[1]]][y]), 
            paste(c, " = ", param.vals[[c[1]]][z]), "\n") 
      }
    } 
  } #close param.vals loops
  
  return(list(out.array, out.qual, 
              temp.params[-which(names(temp.params) %in% c(a, b, c))]))
}
