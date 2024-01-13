## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Model functions
## JM, Last updated 2024-01-13

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

## Functions for determining qualitative behaviors.=============================
# For each parameter set (i.e., row in the parameter matrix), run the ODE and 
#  use the last 10% of timesteps ('last.bit') to characterize qualitative 
#  outcomes for each state variable: positive (TRUE/FALSE), stable equilibrium 
#  (TRUE/FALSE), cyclic dynamics (according to 2 different diagnostics).
# Arguments:
#  param.row = row of parameter dataframe to use;
#  param.df = dataframe with each row specifying a set of parameter values to use;
#  yini = vector containing labeled initial values of state variables;
#  times = vector of timesteps.

run_qual <- function(param.row, param.df, yini, times) {
  temp.params <- param.df[param.row,]
  out <- deSolve::ode(y = yini, times = times, func = full, parms = temp.params)
  out <- data.frame(out)
  max.row <- min( max(out$time), min( which(is.na( out[,labels(yini)] )) ) )
  last.bit <- out[ (0.9 * max.row) : (max.row - 1), ]
  
  var_pos <- sapply(last.bit[,labels(yini)], function(x) { all( x > 1e-08 ) })
  var_neg <- sapply(out[,labels(yini)], function(x) { any( x < 0 ) })
  var_pt.eq <- sapply(last.bit[,labels(yini)], function(x) { isTRUE( max(x) - min(x) < 1e-03 ) })
  var_cyc1 <- sapply(last.bit[,labels(yini)], function(x) { 
    isTRUE( abs(which(x == min(x)) - which(x == max(x))) < length(last.bit[,1]) - 1) })
  var_cyc2 <- sapply(last.bit[,labels(yini)], function(x) { sd(x) / mean(x) })
  var_NA <- sapply(last.bit[,labels(yini)], function(x) { any( is.na(x) == TRUE ) })
  
  return(c(var_pos, var_pt.eq, var_cyc1, var_cyc2, var_neg, var_NA))
}

# Function to run 'run.qual' across combinations of two parameters (defined in 
#  'x.axis' and 'y.axis' arguments) across given values (defined in 'x.range' 
#  and 'y.range' arguments).
# The 'exactY' argument controls whether the function assigns y.range values to 
#  all parameters that contain the text defined in y.axis (this is great for,
#  e.g., using 'ep' to change ep_s and ep_p simultaneously; exactY = FALSE) or 
#  only to parameters that EXACTLY MATCH the text defined in y.axis (e.g., for 
#  using 'm' to change m without also changing mu; exactY = TRUE).
# Other arguments:
#  n.timesteps = number of timesteps;
#  param.vec = vector of parameter values;
#  by = timestep interval;
#  num.cores = number of cores to use for running the function.

qual_out <- function( n.timesteps, param.vec, x.axis, x.range, y.axis, y.range, 
                      by = 1, exactY = FALSE, num.cores = 1 ) {
  
  # Set up a reference dataframe with a row for each combination of parameter
  #  values given by x.range and y.range.
  axes <- expand.grid(x.range, y.range)
  # Set up a parameter dataframe with each row containing all parameter values
  #  defined in 'param.vec', with x.range*y.range rows.
  param.df <- as.data.frame(lapply(param.vec, rep, dim(axes)[1]))
  
  # Replace values of parameters 'x.axis' and 'y.axis' in the parameter dataframe.
  param.df[,which( grepl(x.axis, colnames(param.df)) )] <- axes[,1]
  
  if(exactY == TRUE) { param.df[,y.axis] <- axes[,2] }
  if(exactY == FALSE) { param.df[,which( grepl( y.axis, colnames(param.df) ) )] <- axes[,2] }
  
  times <- seq(0, n.timesteps, by = by) #set timesteps for simulation
  yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001) #set initial conditions
  
  # Set up output dataframe with columns of 'var_pos', 'var_pt.eq', 'var_cyc1',
  #  'var_cyc2', 'var_neg', and 'var_NA' for each of the five state variables.
  var_pos <- sapply( labels(yini), FUN = function(x) { paste0(x, ".pos") } )
  var_pt.eq <- sapply( labels(yini), FUN = function(x) { paste0(x, ".pt.eq") } )
  var_cyc1 <- sapply( labels(yini), FUN = function(x) { paste0(x, ".cyc1") } )
  var_cyc2 <- sapply( labels(yini), FUN = function(x) { paste0(x, ".cyc2") } )
  var_neg <- sapply( labels(yini), FUN = function(x) { paste0(x, ".neg") } )
  var_NA <- sapply( labels(yini), FUN = function(x) { paste0(x, ".NA") } )
  
  out <- data.frame( matrix( nrow = dim(axes)[1], ncol = 30 ) )
  colnames(out) <- c( var_pos, var_pt.eq, var_cyc1, var_cyc2, var_neg, var_NA )
  out <- data.frame( cbind( param.df, out ) ) # combine parameter and output dataframes
  start.col <- which(colnames(out) == "S.pos")
  stop.col <- which(colnames(out) == "H.NA")
  
  # Use 'run_qual' function to populate output dataframe  with qualitative 
  #  behaviors for each parameter set.
  out[1:dim(out)[1],start.col:stop.col] <- t(simplify2array(
    parallel::mclapply(X = 1:dim(out)[1], FUN = run_qual, 
                       param.df = param.df, yini = yini, 
                       times = times, mc.cores = num.cores )))
  
  # Summarize qualitative outcomes: positive point equilibrium (all state variables),
  #  extinction (all state variables), cycles (all state variables), or NA (none of the
  #  above, including differences in qualitative outcomes among state variables).
  out <- data.frame(out)
  out[,c('sum.pos', 'sum.pt.eq', 'sum.cyc1', 'sum.neg', 'sum.NA')] <- sapply( 
    list( out[,var_pos], out[,var_pt.eq], out[,var_cyc1],
          out[,var_neg], out[,var_NA] ), rowSums )
  out$sum.pos.pt.eq <- with(out, sum.pos + sum.pt.eq)
  out$qual <- with(out, 
                   ifelse(sum.NA == 5, "unk",
                          ifelse(sum.pos.pt.eq == 10, "pos.pt.eq", 
                                 ifelse( sum.pos == 0, "extinct",
                                         ifelse(sum.cyc1 == 5 & sum.pos == 5, "cycles", NA )))))
  
  return(out)
}



## Functions for quantifying equilibrium densities.=============================
# For each combination of parameter values, run the ODE across many timesteps and 
#  quantify equilibrium densities for each state variable.
# Requires defined objects 'times' (vector of timesteps), 'params' (vector of
#  named parameter values), and 'param.vals' (list including vectors named a, b, 
#  and c, each containing range of parameter values).

funcy <- function(a, b, c) {
  
  # Set up output array.
  times <- get('times')
  
  out.array <- array( NA, dim = c( length( param.vals[[a[1]]] ),
                                   length( param.vals[[b[1]]] ), 
                                   length( param.vals[[c[1]]] ), 
                                   (0.1*length(times)), 5) )
  all.names <- list( param.vals[[a[1]]], param.vals[[b[1]]], param.vals[[c[1]]], 
                       1:(0.1*length(times)), c('S', 'P', 'I_s', 'I_p', 'H') )
  names(all.names) <- list( a[1], b[1], c[1], 'timestep', 'variable' )
  dimnames(out.array) <- all.names
  
  # For each combination of parameter values, solve ODE and save equilibrium densities.
  for ( x in 1:length(param.vals[[a[1]]]) ) {
    for ( y in 1:length(param.vals[[b[1]]]) ) {
      for ( z in 1:length(param.vals[[c[1]]]) ) {
        
        # Set parameter values using values of a, b, and c.
        temp.params <- get('params')
        temp.params[a] <- param.vals[[a[1]]][x]
        temp.params[b] <- param.vals[[b[1]]][y]
        temp.params[c] <- param.vals[[c[1]]][z]
        
        skip_to_next <- FALSE
        
        # Set initial values for state variables and use ODE solver to get
        #  densities through time.
        yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001)
        out <- deSolve::ode(y = yini, times = times, func = full, parms = temp.params)
        to.length <- min(dim(out)[1], 0.1*length(times))
        
        # Save equilibrium densities to 'out.array' (or skip if error encountered).
        tryCatch({
          out.array[x,y,z,,] <- out[(dim(out)[1]-(to.length-1)):(dim(out)[1]), 
                                    c('S', 'P', 'I_s', 'I_p', 'H')]
        }, error = function(e) { skip_to_next <<- TRUE })
        
        if(skip_to_next) { next }  
        
        # Print parameter values used in current run of loop to track function progress.
        cat(paste(a, " = ", param.vals[[a[1]]][x]), 
            paste(b, " = ", param.vals[[b[1]]][y]), 
            paste(c, " = ", param.vals[[c[1]]][z]), "\n") 
      }
    } 
  } #close param.vals loops
  
  return(reshape2::melt(out.array))
}
