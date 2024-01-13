## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Simulations and figures of cycling behavior
## JM, Last updated 2024-01-04

# Source model functions for model simulations.
source("./inducible-defense-model-v0/01_model-functions.R")

# Load required packages.
require(ggplot2)
require(gridExtra)
require(reshape2)
require(ggtext)

# Set theme for figures.
theme_cyc <- function() { theme_bw() + theme(
  panel.grid=element_blank(),
  axis.text = element_text(size = 14, color = "black"),
  axis.title = element_text(size = 16, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(3,15,3,3),
  panel.border = element_blank(), 
  axis.line = element_line(),
  legend.position = "none"
) }

## Figure 2.====================================================================
## Summary of cycling behavior as a function of theta and high r at low m/low g:
## Panels A-C: bifurcation diagrams of all state variables ~ theta;
## Panels D-E: densities through time at theta = 1.5 & theta = 4;
## Panels F-G: cycle trajectories in the S + P by H phase plane at theta = 1.5 & theta = 4.

## (A) Simulations of cycling behavior.

# Initial parameter list.
params <- c(
  b = 6,
  theta = 1.5,
  g = 0.15,
  x = 1,
  tau_s = 0.8,
  tau_p = 0.8,
  c = 0.8,
  a = 0.1,
  ep_s = 27,
  ep_p = 27,
  w = 1,
  r = 8,
  mu = 0.15,
  m = 0.05
)

times <- seq(0,500, by = 0.5)

# Set list of parameters to vary.
param.vals <- list(
  'theta' = seq(0, 5, by = 0.1),
  'tau_s' = 0.8,
  'tau_p' = 0.8,
  'ep_s' = 27,
  'ep_p' = 27
)

# Run simulation across values specified in 'param.vals'
out1 <- funcy('theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
out2 <- reshape2::dcast(out1, theta + tau_s + ep_s + timestep ~ variable, margin = "value")
# Add columns for prop. of seedlings, mean induction of seedlings, and mean induction
# of mature plant biomass
out2[c('S.prop', 'i_s_mean', 'i_p_mean')] <- c( with(out2, S / (S + P)),
                                                with(out2, I_s / S),
                                                with(out2, I_p / P) )

out.bifur <- reshape2::melt(out2, measure.vars = c("S", "P", "i_s_mean", "i_p_mean", "H"))

# Simulations through time at theta = 1.5 & theta = 4.
yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001)
times <- seq(0, 10000, by = 0.01)

params['theta'] <- 1.5
theta1.5 <- deSolve::ode(y = yini, times = times, func = full, parms = params)
theta1.5 <- data.frame(theta1.5)

params['theta'] <- 4
theta4 <- deSolve::ode(y = yini, times = times, func = full, parms = params)
theta4 <- data.frame(theta4)

# Add column to each dataframe for mean induction across all plant biomass.
theta1.5$i_mean <- with(theta1.5, (I_s + I_p) / (S + P))
theta4$i_mean <- with(theta4, (I_s + I_p) / (S + P))

## (B) Plotting.

# Set up environment for bifurcation plots.
cyc_bifur <- ggplot( ) +
  geom_point(aes( x = theta, y = value, color = variable), size = 0.005) +
  theme_cyc() + 
  labs(x = expression(theta), y = "", title = "<b>A.</b>") +
  scale_x_continuous(expand = c(0,0), limits = c(0,4.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.2), 
                     breaks = seq(0,1.5,by = 0.5)) +
  scale_color_manual(values = c("gray85", "black")) 

# Plot bifurcation diagrams for S and P, I_s and I_p, and H ~ theta (Panels A-C).
cyc_bifur_plots <- grid.arrange(
  cyc_bifur %+% subset(out.bifur, (variable == "S" | variable == "P")) +
    annotate("text", x = 2.7, y = 1.09, label = "P", size = 4.5) + 
    annotate("text", x = 2.7, y = 0.73, label = "S", size = 4.5),
  cyc_bifur %+% subset(out.bifur, (variable == "i_s_mean" | variable == "i_p_mean")) +
    labs(x = expression(theta), y = "", title = "<b>B.</b>" ) +
    scale_y_continuous(limits = c(0, 0.3), expand = c(0,0), breaks = c(0,0.1,0.2,0.3))+
    annotate("text", x = 3, y = 0.195, label = expression(bar(I)[P]), parse = TRUE, size = 4.5) +
    annotate("text", x = 3, y = 0.11, label = expression(bar(I)[S]), parse = TRUE, size = 4.5),
  cyc_bifur %+% subset(out.bifur, (variable == "H")) + 
    labs(x = expression(theta), y = "", title = "<b>C.</b>") +
    annotate("text", x = 3, y = 0.8, label = "H", parse = TRUE, size = 4.5) +
    scale_color_manual(values = c("black")) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)),
  nrow = 1)

# Set up environment for plots of state variables through time.
cyc_time <- ggplot( ) +
  geom_line(aes(x = time, y = S), color = "gray85") +
  geom_line(aes(x = time, y = P), color = "black") +
  geom_line(aes(x = time, y = S+P), color = "gray50") +
  geom_line(aes(x = time, y = H), color = "red3", linetype = "dashed") +
  geom_line(aes(x = time, y = ((I_s + I_p) / (S + P))*4), color = "blue", 
            linetype = "dotdash") +
  theme_cyc() + 
  labs(x = "time", y = "Plant biomass or\n herbivore population size",
       title =  "<b>D.</b> &theta; = 1.5") +
  scale_x_continuous(expand = c(0,0), breaks = seq(1000, 1030, by = 10),
                     limits = c(997, 1037)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2), breaks = seq(0, 2, by = 1),
                     sec.axis = sec_axis(~ . / 4, breaks = c(0,0.25,0.5),
                                         name = "Mean induction"))

# Plot state variables through time for theta = 1.5 and theta = 4 (Panels D-E).
cyc_time_plots <- grid.arrange(
  
  cyc_time %+% subset(theta1.5, time >= 997 & time <= 1030) +
    annotate("text", x = 1033, y = 1.65, label = "S + P", size = 4.5) +
    annotate("text", x = 1033, y = 1.1, label = "P", size = 4.5) +
    annotate("text", x = 1033, y = 0.5, label = "S", size = 4.5) +
    annotate("text", x = 1033, y = 0.8, label = expression(paste(bar(I))), 
             size = 4.5, parse = TRUE) +
    annotate("text", x = 1033, y = 0.23, label = "H", size = 4.5),
  
  cyc_time %+% subset(theta4, time >= 997 & time <= 1031) +
    annotate("text", x = 1033, y = 1.05, label = "S + P", size = 4.5) +
    annotate("text", x = 1033, y = 0.62, label = "P", size = 4.5) +
    annotate("text", x = 1033, y = 0.4, label = "S", size = 4.5) +
    annotate("text", x = 1033, y = 0.85, label = expression(paste(bar(I))), 
             size = 4.5, parse = TRUE) +
    annotate("text", x = 1033, y = 0.26, label = "H", size = 4.5) +
    labs(title =  "<b>E.</b> &theta; = 4") +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), 
                       breaks = seq(0,1,by = 0.5),
                       sec.axis = sec_axis(~ . / 4, 
                                           name = "Mean induction")),
  
  nrow = 1)

# Plot cycle trajectory for theta = 1.5 in the S + P by H phase plane (Panel F).
cyc_traj1 <- ggplot(data = subset(theta1.5, time > 1010 & time < 1025)) +
  geom_path(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),  size = 0.8) +
  geom_path(data = theta1.5[99449:99450,], aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),
            size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  geom_path(data = theta1.5[99999:100000,], aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),
            size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  theme_cyc() + theme(
    legend.title = element_text(size = 14, family = "serif"),
    legend.text = element_text(size = 12),
    legend.justification = c(1,0.5),
    legend.position = c(1,0.6),
    legend.key.width = unit(0.95, "lines"), legend.key.height = unit(1, "lines")
  ) +
  geom_point( data = theta1.5[98000 + which(theta1.5[98000:99400,'i_mean'] == min(theta1.5[98000:99400,'i_mean'])),],
              aes(x = S + P, y = H, color = (I_s + I_p)/(S + P)), size = 1.8)+
  geom_point( data = theta1.5[98000 + which(theta1.5[98000:99400,'i_mean'] == max(theta1.5[98000:99400,'i_mean'])),],
              aes(x = S + P, y = H, color = (I_s + I_p)/(S + P)), size = 1.8)+
  scale_x_continuous(limits = c(1,2.4), expand = c(0,0), breaks = c(1,2)) +
  scale_y_continuous(limits = c(0,1.3), expand = c(0,0), breaks = c(0,0.5,1)) +
  labs(x = "S + P", y = "H", color = expression(bar(I)), title = "<b>F.</b> &theta; = 1.5" )+
  scale_color_gradient(low = "gray90", high = "black", limits = c(0.16,0.245), breaks = seq(0.16,0.24, by = 0.04))

# Plot cycle trajectory for theta = 4 in the S + P by H phase plane (Panel G).
cyc_traj2 <- ggplot(data = subset(theta4, time > 1005 & time < 1012)) +
  geom_path(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),  size = 0.8)+
  geom_path(data = theta4[99349:99350,],
            aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)), size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  geom_path(data = theta4[99549:99550,],
            aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)), size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  theme_cyc() + 
  geom_point( data = theta4[98000 + which(theta4[98000:99400,'i_mean'] == min(theta4[98000:99400,'i_mean'])),],
              aes(x = S + P, y = H, color = (I_s + I_p)/(S + P)), size = 1.8)+
  geom_point( data = theta4[98000 + which(theta4[98000:99400,'i_mean'] == max(theta4[98000:99400,'i_mean'])),],
              aes(x = S + P, y = H, color = (I_s + I_p)/(S + P)), size = 1.8) +
  scale_x_continuous(limits = c(0.9,1.1), expand = c(0,0),breaks=seq(0.9,1.1, by = 0.1)) +
  scale_y_continuous(limits = c(0.2,0.5), expand = c(0,0), breaks = seq(0.2, 0.5, by = 0.1)) +
  labs(x = "S + P", y = "H", color = expression(bar(I)), title = "<b>G.</b> &theta; = 4" )+
  scale_color_gradient(low = "gray90", high = "black", limits = c(0.16,0.245))

cyc_traj_plots <- grid.arrange(cyc_traj1, cyc_traj2, nrow = 1)

Fig2 <- grid.arrange(cyc_bifur_plots, cyc_time_plots, cyc_traj_plots,
                     nrow = 3, heights = c(1,1.1,1))

## Figure 3.====================================================================
## Bifurcation of S, P, and H as a function of seedling maturation rate (g) for 
##  increased theta at low m/low g and high m/low g.

## (A) Simulations of cycling behavior.

# Initial parameter list.
params <- c(
  b = 6,
  theta = 4,
  g = 0.15,
  x = 1,
  tau_s = 0.8,
  tau_p = 0.8,
  c = 0.8,
  a = 0.1,
  ep_s = 27,
  ep_p = 27,
  w = 1,
  r = 2,
  mu = 0.15,
  m = 0.05
)

# Set list of parameters to vary.
param.vals <- list(
  'g' = seq(0.025, 1, by = 0.025),
  'tau_s' = 0.8,
  'tau_p' = 0.8,
  'ep_s' = 27,
  'ep_p' = 27
)

times <- seq(0,500, by = 0.5)

# Simulations for low m/low g.
out1 <- funcy('g', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
out2 <- reshape2::dcast(out1, g + tau_s + ep_s + timestep ~ variable, margin = "value")
out.bifur_low.m_g <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))

# Simulations for high m/low g.
params[c('b', 'm', 'g')] <- c(11, 0.8, 0.15)

out1 <- funcy('g', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
out2 <- reshape2::dcast(out1, g + tau_s + ep_s + timestep ~ variable, margin = "value")
out.bifur_high.m_g <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))

## (B) Plotting.

# Set up environment for bifurcation plots.
cyc1 <- ggplot() +
  geom_point(size = 0.03, position = position_dodge(width = 0.01)) +
  stat_summary(fun = "mean", geom = "line") +
  theme_cyc() + 
  labs(y = "") +
  scale_color_manual(values = c("gray85", "black")) 

Fig3 <- grid.arrange(
  # For m = 0.05, g = 0.15
  cyc1 %+% subset(out.bifur_low.m_g, (variable == "S" | variable == "P")) + 
    aes(x = g, y = value, color = variable) +
    geom_vline(xintercept = 0.15, linetype = "dashed") +
    annotate("text", x = 0.85, y = 1, label = "P", size = 4.5) +
    annotate("text", x = 0.85, y = 0.3, label = "S", size = 4.5) +
    labs(x = "g", y = "", title = "<b>A.</b> <i>m</i> = 0.05, <i>b</i> = 6") +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.8), 
                       breaks = seq(0,1.5,by = 0.5)) ,
  
  cyc1 %+% subset(out.bifur_low.m_g, (variable == "H")) + 
    aes(x = g, y = value, color = variable) +
    geom_vline(xintercept = 0.15, linetype = "dashed") +
    annotate("text", x = 0.85, y = 0.55, label = "H", size = 4.5) +
    labs(x = "g", y = "", title = "<b>B.</b> <i>m</i> = 0.05, <i>b</i> = 6") +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.8), 
                       breaks = seq(0,1.5,by = 0.5)) +
    scale_color_manual(values = c("black")) ,
  
  # For m = 0.8, g = 0.15
  cyc1 %+% subset(out.bifur_high.m_g, (variable == "S" | variable == "P")) + 
    aes(x= g, y = value, color = variable) +
    geom_vline(xintercept = 0.15, linetype = "dashed") +
    annotate("text", x = 0.52, y = 0.15, label = "P", size = 4.5) +
    annotate("text", x = 0.52, y = 0.45, label = "S", size = 4.5) +
    labs(x = "g", y = "", title = "<b>C.</b> <i>m</i> = 0.80, <i>b</i> = 11") +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), 
                       breaks = seq(0,1,by = 0.5)) ,
  
  cyc1 %+% subset(out.bifur_high.m_g, (variable == "H")) + 
    aes(x = g, y = value, color = variable) +
    geom_vline(xintercept = 0.15, linetype = "dashed") +
    annotate("text", x = 0.52, y = 0.5, label = "H", size = 4.5) +
    labs(x = "g", y = "", title = "<b>D.</b> <i>m</i> = 0.80, <i>b</i> = 11") +
    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 1.1), 
                       breaks = seq(0,1,by = 0.5)) +
    scale_color_manual(values = c("black")) ,
  
  nrow = 2)


## Figure S1.===================================================================
## Bifurcation diagram of S + P ~ theta at high r for (A-C) increasing x with 
##  constant g and (D-F) increasing x with constant gx.

## (A) Simulations of cycling behavior.

# Initial parameter list.
params <- c(
  b = 6,
  theta = 1.5,
  g = 0.15,
  x = 1,
  tau_s = 0.8,
  tau_p = 0.8,
  c = 0.8,
  a = 0.1,
  ep_s = 27,
  ep_p = 27,
  w = 1,
  r = 8,
  mu = 0.15,
  m = 0.05
)

times <- seq(0,500, by = 0.5)

# Set list of parameter to vary.
param.vals <- list(
  'theta' = seq(0,5, by = 0.1),
  'tau_s' = 0.8,
  'tau_p' = 0.8,
  'ep_s' = 27,
  'ep_p' = 27
)

# List with combinations of b, x, and g to use.
list1 <- list(
  c(3, 2, 0.15),
  c(2, 3, 0.15),
  c(1, 6, 0.15),
  c(5.55, 2, 0.075),
  c(5.4, 3, 0.05),
  c(5.25, 6, 0.025)
)

# For each element in list1, reassign b, x, and g, find equilibrium behavior as
#  a function of theta, and save output to an object.

for (i in 1:length(list1)) {
  params[c('b', 'x', 'g')] <- list1[[i]]
  out1 <- funcy('theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
  out2 <- reshape2::dcast(out1, theta + tau_s + ep_s + timestep ~ variable, margin = "value")
  out.bifur <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))
  assign(paste0('x_b', params['b'], '_x', params['x']), out.bifur)
}

## (B) Plotting.

# Set up environment for bifurcation plots.
cycx1 <- ggplot( ) +
  geom_point(aes( x = theta, y = value, color = variable), size = 0.015) +
  theme_cyc() +
  labs(x = expression(theta), y = "", title = "<b>A.</b>") +
  scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.3), 
                     breaks = seq(0,1.5,by = 0.5)) +
  scale_color_manual(values = c("gray85", "black")) 

FigS2 <- grid.arrange(
  cycx1 %+% subset(x_b3_x2, (variable == "S" | variable == "P")) +
    labs(title = "<b>A.</b> <i>x</i> = 2, <i>g</i> = 0.075"),
  cycx1 %+% subset(x_b2_x3, (variable == "S" | variable == "P")) +
    labs(title = "<b>B.</b> <i>x</i> = 3, <i>g</i> = 0.05"),
  cycx1 %+% subset(x_b1_x6, (variable == "S" | variable == "P")) +
    labs(title = "<b>C.</b> <i>x</i> = 6, <i>g</i> = 0.025"),
  cycx1 %+% subset(x_b5.55_x2, (variable == "S" | variable == "P")) +
    labs(title = "<b>D.</b> <i>x</i> = 2, <i>g</i> = 0.075"),
  cycx1 %+% subset(x_b5.4_x3, (variable == "S" | variable == "P")) +
    labs(title = "<b>E.</b> <i>x</i> = 3, <i>g</i> = 0.05"),
  cycx1 %+% subset(x_b5.25_x6, (variable == "S" | variable == "P")) +
    labs(title = "<b>F.</b> <i>x</i> = 6, <i>g</i> = 0.025"),
  nrow = 2)


## Figures S3-5.================================================================
## Bifurcation diagrams of S, P, and H as function of herbivore intrinsic growth 
##  rate (r) (Fig. S3), induction sensitivity to damage (c) (Fig. S4), and 
##  induction decay rate (a) (Fig. S5) for increased theta at low m/low g, 
##  high m/low g, and low m/moderate g.

## (A) Simulations of cycling behavior.

# Initial parameter list.
params <- c(
  b = 6,
  theta = 4,
  g = 0.15,
  x = 1,
  tau_s = 0.8,
  tau_p = 0.8,
  c = 0.8,
  a = 0.1,
  ep_s = 27,
  ep_p = 27,
  w = 1,
  r = 2,
  mu = 0.15,
  m = 0.05
)

times <- seq(0,500, by = 0.5)

# Set list of parameters to vary.
param.vals <- list(
  'r' = seq(1,12, by = 0.25),
  'c' = seq(0.2, 5, by = 0.1),
  'a' = seq(0.02, 0.4, by = 0.01),
  'tau_s' = 0.8,
  'tau_p' = 0.8,
  'ep_s' = 27,
  'ep_p' = 27
)

# List with combinations of b, m, and g to use.
list1 <- list(
  c(6, 0.05, 0.15),
  c(11, 0.8, 0.15),
  c(2.43, 0.05, 0.5)
)

# For each element in list1, reassign b, m, and g, find equilibrium behavior
#  as a function of r, c, or a, and save output to an object.

for (i in 1:length(list1)) {
  
  params[c('b', 'm', 'g')] <- list1[[i]]

  out1 <- funcy('r', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
  out2 <- reshape2::dcast(out1, r + tau_s + ep_s + timestep ~ variable, margin = "value")
  out.bifur <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))
  assign(paste0('r_m', params['m'], '_g', params['g']), out.bifur)
  
  out1 <- funcy('c', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
  out2 <- reshape2::dcast(out1, c + tau_s + ep_s + timestep ~ variable, margin = "value")
  out.bifur <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))
  assign(paste0('c_m', params['m'], '_g', params['g']), out.bifur)
  
  out1 <- funcy('a', c('tau_s', 'tau_p'), c('ep_s', 'ep_p'))
  out2 <- reshape2::dcast(out1, a + tau_s + ep_s + timestep ~ variable, margin = "value")
  out.bifur <- reshape2::melt(out2, measure.vars = c("S", "P", "H"))
  assign(paste0('a_m', params['m'], '_g', params['g']), out.bifur)

}

## (B) Plotting.

# (uses 'cyc1' plotting environment defined for Fig. 3 above)

FigS3 <- grid.arrange(
  
  # For m = 0.05, g = 0.15
  cyc1 %+% subset(r_m0.05_g0.15, variable == "S" | variable == "P")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>A.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    annotate("text", x = 10.5, y = 0.75, label = "P", size = 4.5) +
    annotate("text", x = 10.5, y = 0.25, label = "S", size = 4.5) +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3),
  cyc1 %+% subset(r_m0.05_g0.15, variable == "H")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>B.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    annotate("text", x = 10.5, y = 0.45, label = "H", size = 4.5) +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.8, g = 0.15
  cyc1 %+% subset(r_m0.8_g0.15, variable == "S" | variable == "P")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>C.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5), n.breaks = 3),
  cyc1 %+% subset(r_m0.8_g0.15, variable == "H")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>D.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.5), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.05, g = 0.5
  cyc1 %+% subset(r_m0.05_g0.5, variable == "S" | variable == "P")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>E.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3),
  cyc1 %+% subset(r_m0.05_g0.5, variable == "H")+
    aes(x = r, y= value, color = variable) +
    labs(title = "<b>F.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,12)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  nrow = 3)


FigS4 <- grid.arrange(
  
  # For m = 0.05, g = 0.15
  cyc1  %+% subset(c_m0.05_g0.15, variable == "S" | variable == "P") +
    annotate("text", x = 4, y = 0.75, label = "S", size = 4.5) +
    annotate("text", x = 4, y = 0.4, label = "P", size = 4.5) +
    aes(x = c, y= value, color = variable) +
    labs(title = "<b>A.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3),
  cyc1 %+% subset(c_m0.05_g0.15, variable == "H")+
    aes(x = c, y= value, color = variable) +
    annotate("text", x = 4, y = 1.65, label = "H", size = 4.5) +
    labs(title = "<b>B.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.6), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.8, g = 0.15
  cyc1 %+% subset(c_m0.8_g0.15, variable == "S" | variable == "P")+
    aes(x = c, y= value, color = variable) +
    labs(title = "<b>C.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.6), n.breaks = 3),
  cyc1 %+% subset(c_m0.8_g0.15, variable == "H")+
    aes(x = c, y= value, color = variable) +
    labs(title = "<b>D.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2.6), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.05, g = 0.5
  cyc1 %+% subset(c_m0.05_g0.5, variable == "S" | variable == "P")+
    aes(x = c, y= value, color = variable)+
    labs(title = "<b>E.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3),
  cyc1 %+% subset(c_m0.05_g0.5, variable == "H")+
    aes(x = c, y= value, color = variable) +
    labs(title = "<b>F.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,5)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  nrow = 3)

FigS5 <- grid.arrange(
  
  # For m = 0.05, g = 0.15
  cyc1  %+% subset(a_m0.05_g0.15, variable == "S" | variable == "P") +
    aes(x = a, y= value, color = variable) +
    annotate("text", x = 0.2, y = 0.85, label = "P", size = 4.5) +
    annotate("text", x = 0.2, y = 0.45, label = "S", size = 4.5) +
    labs(title = "<b>A.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1.05), n.breaks = 3),
  cyc1 %+% subset(a_m0.05_g0.15, variable == "H")+
    aes(x = a, y= value, color = variable) +
    annotate("text", x = 0.2, y = 1.05, label = "H", size = 4.5) +
    labs(title = "<b>B.</b> <i>m</i> = 0.05, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,2), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.8, g = 0.15
  cyc1 %+% subset(a_m0.8_g0.15, variable == "S" | variable == "P")+
    aes(x = a, y= value, color = variable) +
    labs(title = "<b>C.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,5.5), n.breaks = 3),
  cyc1 %+% subset(a_m0.8_g0.15, variable == "H")+
    aes(x = a, y= value, color = variable) +
    labs(title = "<b>D.</b> <i>m</i> = 0.80, <i>g</i> = 0.15", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,5), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  
  # For m = 0.05, g = 0.5
  cyc1 %+% subset(a_m0.05_g0.5, variable == "S" | variable == "P")+
    aes(x = a, y= value, color = variable)+
    labs(title = "<b>E.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3),
  cyc1 %+% subset(a_m0.05_g0.5, variable == "H")+
    aes(x = a, y= value, color = variable) +
    labs(title = "<b>F.</b> <i>m</i> = 0.05, <i>g</i> = 0.50", y = "") +
    scale_x_continuous(expand = c(0,0), limits = c(0,0.4)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), n.breaks = 3) +
    scale_color_manual(values = c("black")),
  nrow = 3)
