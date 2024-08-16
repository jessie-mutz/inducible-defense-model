## Mutz, J. & K.C. Abbott. Life history modulates effects of inducible defenses 
##  on consumer-resource dynamics.
## Supplemental code: Figures
## JM, Last updated 2024-06-05

# Source model functions.
source("./inducible-defense-model/01_model-functions.R")

# Load packages.
require(ggplot2)
require(gridExtra)
require(ggtext)
require(doBy)
require(reshape2)

## Figure 1.====================================================================
## Panel A: Behavior ~ seedling maturation rate (g) x reproductive plant mortality (m)
##  across values of theta (inducible defense costs; rows) and tau (inducible
##  defense responsiveness; columns).
## Panel B: Threshold value of epsilon needed for population persistence ~ seedling
##  maturation rate (g) x reproductive plant mortality (m).
## Panel C: Effect of increasing seedling maturation rate (g) on model behavior 
##  (stabilizing, destabilizing, or no change) ~ seedling maturation rate (g) x
##  reproductive plant mortality (m).

# Set theme for figures.
theme_qual <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(),
  strip.background = element_blank(),
  axis.text = element_text(size = 11),
  plot.title = ggtext::element_markdown(size = 12)
)
}

## Read in simulation output for Panel A.
qual_df <- read.csv("04_tau-theta-ep27-qual.csv")
qual_df$qual2 <- ifelse(is.na(qual_df$qual) == TRUE, "extinct", qual_df$qual)

# Plot qualitative behavior ~ g x m, across values of tau and theta.
tau_fig <- ggplot(subset(qual_df, (ep_s == 27) & 
                           (theta == 3 | theta == 3.5 | theta == 4)),
                  aes(x = g, y = m, color = qual2, shape = qual2)) +
  geom_hline(yintercept = -Inf, size = 1) +
  geom_vline(xintercept = -Inf, size = 1) +
  geom_point(size = 1.35) + 
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  theme_qual() +
  theme(
    strip.text = element_text(hjust = 0, size = 11),
    axis.title = element_text(size = 12),
    legend.key.height = unit(2, "lines"),
    legend.text = element_text(size = 11.5),
    panel.spacing = unit(1, "lines")
  )  +
  scale_x_continuous(limits = c(0,1), expand = c(0,0), breaks = seq(0,1,by = 0.5),
                     labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = seq(0,1,by = 0.5),
                     labels = c("0", "0.5", "1")) +
  facet_grid(rows = vars(theta), cols = vars(tau_s), 
             labeller = label_bquote(cols = tau == .(tau_s), rows = theta == .(theta))) +
  labs(x = "Seedling maturation rate (g)", y = "Reproductive plant mortality rate (m)",
       color = "", shape = "") +
  scale_color_manual(values = c('paleturquoise', 'firebrick4', 'gray60'),
                     breaks = c('pos.pt.eq', 'cycles', 'extinct'),
                     labels = c('stable', 'cycles', 'extinct')) +
  scale_shape_manual(values = c(17, 15, 16),
                     breaks = c('pos.pt.eq', 'cycles', 'extinct'),
                     labels = c('stable', 'cycles', 'extinct')) 

## Read in simulation output for Panel B.
qual_ep_df <- read.csv("05_tau0.8-theta-ep-qual.csv")
qual_ep_df$qual2 <- ifelse(is.na(qual_ep_df$qual) == TRUE, "extinct", qual_ep_df$qual)

# Find minimum value of epsilon at which population persists for each parameter
#  combination.
qual_ep_df2 <- summaryBy(ep_s ~ rep + g + m + tau_s + theta,
                         FUN = min, 
                         data = subset(qual_ep_df, qual == "pos.pt.eq" |
                                         qual == "cycles"))

# Plot minimum value of epsilon at which population persists ~ g x m.
ep.fig <- ggplot(subset(qual_ep_df2, (theta == 4)),
                 aes(x = g, y = m, color = ep_s.min)) +
  geom_point() + 
  theme_qual() +
  theme(
    plot.margin = margin(5,15,5,5),
    strip.text = element_text(hjust = 0),
    legend.title = ggtext::element_markdown(size = 8),
    legend.text = element_text(size = 8),
    legend.margin = margin(0,0,0,0)
  ) +
  labs(x = "Seedling maturation\nrate (g)",
       y = "Reproductive plant\nmortality rate (m)",
       color = "Induced defense<br>effectiveness (&epsilon;)", 
       tag = "A.") +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_color_gradient(low = "white", high = "black", limits = c(20,40))

ep.fig %+% subset(qual_ep_df2, theta == 2) +
  scale_color_gradient(low = "white", high = "paleturquoise4", limits = c(10,20),
                       breaks = seq(10, 20, by = 2))

## Read in simulation output for Panel C.
qual_g_df <- read.csv("06_tau0.8-theta4-ep-g0.955-qual.csv")
qual_g_df$qual2 <- ifelse(is.na(qual_g_df$qual) == TRUE, "extinct", qual_g_df$qual)
qual_g_df[which(qual_g_df$rep == 46), 'qual2'] <- "cycles"

# Add values of seedling maturation rate (g) for which there are cycles
sub <- subset(qual_ep_df, ep_s == 27 & theta == 4 & qual == "cycles")
qual_g_df$g2 <- sub$g

# Plot qualitative behavior when seedling maturation rate (g) is increased.
g.fig <- ggplot(data = qual_g_df, 
                aes(x = g2, y = m, color = qual2, shape = qual2)) +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_qual() +
  theme(
    plot.margin = margin(5,15,5,5),
    strip.text = element_text(hjust = 0),
    legend.position = c(1,1),
    legend.justification = c(1,1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  labs(x = "Seedling maturation\nrate (g)",
       y = "Reproductive plant\nmortality rate (m)",
       color = "Effect of increased g\non equilibrium behavior",
       shape = "Effect of increased g\non equilibrium behavior") +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_color_manual(values = c('paleturquoise', 'firebrick4', 'gray60'),
                     breaks = c("pos.pt.eq", "cycles", "extinct"),
                     labels = c("Stabilizing\n(cycles to stable)",
                                "None",
                                "Destabilizing\n(cycles to extinct)")) +
  scale_shape_manual(values = c(16, 15, 17),
                     breaks = c("pos.pt.eq", "cycles", "extinct"),
                     labels = c("Stabilizing\n(cycles to stable)",
                                "None",
                                "Destabilizing\n(cycles to extinct)"))

## Combine panels.
Fig1 <- grid.arrange(
  tau_fig +
    labs(color = "", shape = "", tag = "A."),
  grid.arrange(
    ep.fig +   
      labs(color = "Threshold<br>inducible defense<br>effectiveness (&epsilon;)", 
           tag = "B.", title = "&theta; = 4"),
    g.fig +  
      labs(color = "Effect of increased g\non equilibrium behavior",
           shape = "Effect of increased g\non equilibrium behavior",
           title = "&theta; = 4", tag = "C."),
    nrow = 1, widths = c(1, 0.8)
  ),
  nrow = 2, heights = c(1, 0.75))


## Figure 2.====================================================================
## Panel A: Behavior ~ seedling maturation rate (g) x reproductive plant mortality (m)
##  for herbivore intrinsic growth rate (r) = 8, across values of theta (inducible
##  defense costs)
## Panel B: Plant biomass or herbivore population size through time, across values
##  of theta (inducible defense costs)
## Panel C: Cycle trajectories in the S + P by H phase plane, across values of theta
##  (inducible defense costs)

# Set theme for figures.
theme_cyc <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  panel.spacing = unit(1, "lines"),
  panel.border = element_blank(), 
  axis.line = element_line(),
  axis.text = element_text(size = 11, color = "black"),
  axis.title = element_text(size = 12),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(10,10,10,0),
  legend.position = "top",
  legend.margin = margin(0,0,0,0),
  legend.background = element_blank(),
  legend.location = "plot",
  strip.background = element_blank(),
  strip.text = element_blank()
) }

## Read in simulation output for Panel A.
qual_r8_df <- read.csv("07_tau-theta-ep27-r8-qual.csv")
qual_r8_df$qual2 <- ifelse(is.na(qual_r8_df$qual) == TRUE, "extinct", qual_r8_df$qual)

# Plot qualitative behavior ~ g x m at high inducible defense costs.
r8_1 <- ggplot(subset(qual_r8_df, (ep_s == 27) & tau_s == 0.8 &
                        (theta == 2 | theta == 3 | theta == 4)),
               aes(x = g, y = m, color = qual2, shape = qual2)) +
  geom_hline(yintercept = -Inf, size = 1) +
  geom_vline(xintercept = -Inf, size = 1) +
  geom_point(size = 1.2) + 
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  theme_cyc() +
  theme(
    legend.text = element_text(size = 11),
    legend.justification.top = "left",
    legend.key = element_blank()
  ) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0), breaks = seq(0,1,by = 0.5),
                     labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), breaks = seq(0,1,by = 0.5),
                     labels = c("0", "0.5", "1")) +
  facet_wrap(~theta, labeller = label_bquote(theta == .(theta)), nrow = 3) +
  labs(x = "Seedling maturation\nrate (g)", 
       y = "Reproductive plant mortality rate (m)",
       color = "", shape = "", tag = "A.") +
  scale_color_manual(#values = c('lightblue1', 'tomato4', 'gray75'),
    values = c('paleturquoise', 'firebrick4', 'gray60'),
    breaks = c('pos.pt.eq', 'cycles', 'extinct'),
    labels = c('stable', 'cyc.', 'ext.')) +
  scale_shape_manual(values = c(17, 15, 16),
                     breaks = c('pos.pt.eq', 'cycles', 'extinct'),
                     labels = c('stable', 'cyc.', 'ext.'))

## Use one set of seedling maturation rate (g) and reproductive plant mortality 
##  rate (m) to simulate initial dynamics.

# Initial parameter list.
params <- c(
  b = qual_r8_df[which(qual_r8_df$rep == 402),'b'][1],
  theta = 1.5,
  g = qual_r8_df[which(qual_r8_df$rep == 402),'g'][1],
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
  m = qual_r8_df[which(qual_r8_df$rep == 402),'m'][1]
)

yini <- c(S = 0.05, P = 0.05, I_s = 0, I_p = 0, H = 0.001)
times <- seq(0, 10000, by = 0.01)

# Simulate initial dynamics for theta = 2.
params['theta'] <- 2
theta2 <- deSolve::ode(y = yini, times = times, func = full, parms = params)
theta2 <- data.frame(theta2)

# Simulate initial dynamics for theta = 3.
params['theta'] <- 3
theta3 <- deSolve::ode(y = yini, times = times, func = full, parms = params)
theta3 <- data.frame(theta3)

# Simulate initial dynamics for theta = 4.
params['theta'] <- 4
theta4 <- deSolve::ode(y = yini, times = times, func = full, parms = params)
theta4 <- data.frame(theta4)

## For Panel B: Select simulation output for timesteps from 0-100.
theta_out <- rbind(
  cbind(subset(theta2, time < 100), theta = 2),
  cbind(subset(theta3, time < 100), theta = 3),
  cbind(subset(theta4, time < 100), theta = 4)
)

# Add columns for mean inducible defense level across both seedling and 
#  reproductive plant biomass; multiply by 4 for plotting.
theta_out$I_mean <- with(theta_out, (I_s + I_p) / (S + P))
theta_out$I_mean_sc <- with(theta_out, I_mean*4)
theta_out$SP_tot <- with(theta_out, S + P)

# Melt dataframe so that variables to be plotted are all within a single column.
theta_out_melt <- melt(data = theta_out, 
                       measure.vars = c('S', 'P', 'I_mean_sc', 'H', 'SP_tot'),
                       id.vars = c('time', 'theta'))

# Plot abundances through time for each value of inducible defense costs.
cyc_time <- ggplot( theta_out_melt) +
  geom_hline(yintercept = -Inf) +
  geom_line(aes(x = time, y = value, color = variable,
                linetype = variable), size = 0.6) +
  theme_cyc() +   
  theme(
    legend.text = element_text(size = 11.5),
    legend.justification.top = "left",
    legend.key = element_blank()
  ) +
  labs(x = "time\n ", y = "Plant biomass or herbivore pop. size",
       color = "", linetype = "", tag = "B.") +
  scale_x_continuous(expand = c(0,0), limits = c(0,80)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 3), breaks = seq(0, 3, by = 1),
                     sec.axis = sec_axis(~ . / 4, breaks = c(0, 0.2, 0.4, 0.6),
                                         name = "Mean induction")) +
  scale_color_manual(values = c('gray85', 'black', 'gray50', 'slateblue', 'chartreuse2'),
                     breaks = c('S', 'P', 'SP_tot', 'H', 'I_mean_sc'),
                     labels = c('S', 'P', 'S + P', 'H', expression(bar(I)))) +
  scale_linetype_manual(values = c('solid', 'solid', 'solid', 'dashed', 'dotdash'),
                        breaks = c('S', 'P', 'SP_tot', 'H', 'I_mean_sc'),
                        labels = c('S', 'P', 'S + P', 'H', expression(bar(I)))) +
  facet_wrap(~ theta, labeller = label_bquote(theta == .(theta)), ncol = 1) 

## For Panel C: Select simulation output across the length of one cycle.
theta_out <- rbind(
  cbind(subset(theta2, time > 1004 & time < 1016.2), theta = 2),
  cbind(subset(theta3, time == 100), theta = 3),
  cbind(subset(theta4, time > 1005 & time < 1012), theta = 4)
)

# Plot trajectories at equilibrium within the plant biomass x herbivore 
#  population size phase plane. 
traj1 <- ggplot(data = theta_out) +
  geom_hline(yintercept = -Inf) +
  geom_path(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),  size = 0.6) +
  geom_path(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),
            size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  geom_path(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P)),
            size = 0.6, arrow = arrow(length = unit(1, "lines"), angle = 20)) +
  geom_point(aes(x = S + P, y = H, color = (I_s + I_p) / (S + P), 
                 size = as.factor(theta))) +
  theme_cyc() +
  theme(    
    legend.text = element_text(size = 10.5),
    legend.key.width = unit(2, "lines"),
    legend.key.height = unit(0.5, "lines"),
    legend.title = element_text(size = 12.5, vjust = 1)
  ) +
  labs(x = "Total plant biomass\n(S + P)", y = "Herbivore pop. size (H)",
       color = expression(bar(I)), tag = "C.") +
  facet_wrap(~ theta, labeller = label_bquote(theta == .(theta)),
             ncol = 1) +
  scale_x_continuous(limits = c(0.8, 1.8), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  scale_size_manual(values = c(0, 2.5, 0), guide = "none") +
  scale_color_gradient(low = "gray95", high = "black", 
                       limits = c(0.17, 0.23), breaks = seq(0.18, 0.22, by = 0.02))

## Combine panels.
Fig2 <- grid.arrange(
  r8_1 + 
    annotate("point", x = 0.149, y = 0.067, shape = 21, 
             color = "gold", size = 3.5, stroke = 1.3), 
  cyc_time, 
  traj1, 
  nrow = 1, widths = c(1, 1.3, 1))


## Figure 3.====================================================================

# Set theme for figures.
theme_quan <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 14, color = "black"),
  axis.title.y = element_text(size = 16, family = "serif"),
  axis.title.x = ggtext::element_markdown(size = 16, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(3,15,3,3),
  panel.border = element_blank(), 
  axis.line = element_line(),
  legend.position = "none"
) }

# Read in simulation output.
quan_df <- read.csv("08_tau-theta-ep-quan.csv")

# Designate simulation runs where the system went to extinction.
quan_df[,'test'] <- apply(quan_df[,15:24], MARGIN = 1, 
                          function(x) { any(x < 0)})

# Select only four combinations of m and g to plot, and remove any extinctions.
quan_df <- subset(quan_df, timestep == 100 & test == FALSE &
                    (m == 0.05 | m == 0.8) & (g == 0.15 | g == 0.5))

# Melt dataframe to have variables in single column.
quan_df_melt <- melt(data = quan_df, 
                     measure = c('S', 'P', 'I_s', 'I_p', 'H', 
                                 'i_s_mean', 'i_p_mean',
                                 'P.prop', 'H_SP', 'seeds'))

rm(quan_df) # remove unmelted df -- very big and no longer needed

# Set up environment for plotting equilibrium as a function of tau, with lines
#  for each combination of m and g values.
quan_fig1 <- ggplot( ) +
  geom_line(aes( x = tau_s, y = value, color = as.factor(g), linetype = as.factor(m) ))+
  theme_quan() + 
  labs(x = "", 
       y = expression(paste( hat(S) )), 
       title = "<b>A.</b> Seedling biomass")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_color_manual(values = c("gray", "black")) +
  scale_linetype_manual(values = c("solid", "dashed")) 

# Plotting for baseline theta.
quan_base <- grid.arrange( 
  quan_fig1 %+% subset(quan_df_melt, variable == "S" & theta == 2 & ep_s == 30) +
    annotate("segment", x = 0.5, xend = 0.6, y = 4, yend = 4) +
    annotate("segment", x = 0.5, xend = 0.6, y = 3.5, yend = 3.5, linetype = "dashed") +
    annotate("point", x = 0.55, y = 3, shape = 22, fill = "gray80", color = "white", size = 5.5) +
    annotate("point", x = 0.55, y = 2.5, shape = 22, fill = "black", color = "white", size = 5.5) +  
    annotate("text", x = 0.63, y = 4, label = expression(paste(italic(m), " = 0.05")), hjust = 0) +
    annotate("text", x = 0.63, y = 3.5, label = expression(paste(italic(m), " = 0.80")), hjust = 0) +
    annotate("text", x = 0.63, y = 3, label = expression(paste(italic(g), " = 0.15")), hjust = 0) +
    annotate("text", x = 0.63, y = 2.5, label = expression(paste(italic(g), " = 0.50")), hjust = 0),
  quan_fig1 %+% subset(quan_df_melt, variable == "P" & theta == 2 & ep_s == 30) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(y = expression(paste( hat(P) )), 
         title = "<b>B.</b> Mature plant biomass"),
  quan_fig1 %+% subset(quan_df_melt, variable == "H" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste( hat(H) )), 
         title = "<b>C.</b> Herbivore pop. size" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  quan_fig1 %+% subset(quan_df_melt, variable == "seeds" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste("b", hat(P), "(1 - ", theta, hat(bar(I))[P],")")), 
         title = "<b>D.</b> Seed production" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,7)),
  nrow = 1)

# Plotting for increased theta.
quan_theta <- grid.arrange( 
  quan_fig1 %+% subset(quan_df_melt, variable == "S" & theta == 4 & ep_s == 30) +
    labs(x = "Induction<br>responsiveness&nbsp;(&tau;)",
         y = expression(paste( hat(S) )), 
         title = "<b>I.</b> Seedling biomass"), 
  quan_fig1 %+% subset(quan_df_melt, variable == "P" & theta == 4 & ep_s == 30) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = "Induction<br>responsiveness&nbsp;(&tau;)",
         y = expression(paste( hat(P) )), 
         title = "<b>J.</b> Mature plant biomass"),
  quan_fig1 %+% subset(quan_df_melt, variable == "H" & theta == 4 & ep_s == 30) +
    labs(x = "Induction<br>responsiveness&nbsp;(&tau;)",
         y = expression(paste( hat(H) )), 
         title = "<b>K.</b> Herbivore pop. size" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  quan_fig1 %+% subset(quan_df_melt, variable == "seeds" & theta == 4 & ep_s == 30) +
    labs(x = "Induction<br>responsiveness&nbsp;(&tau;)",
         y = expression(paste("b", hat(P), "(1 - ", theta, hat(bar(I))[P],")")), 
         title = "<b>L.</b> Seed production" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,7)),
  nrow = 1)

# Plotting for variables with no effect of theta.
quan_both <- grid.arrange(
  quan_fig1 %+% subset(quan_df_melt, variable == "i_s_mean" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste( hat(bar(I))[S] )), 
         title = "<b>E.</b> Seedling induction" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.25), breaks = seq(0,0.2,by = 0.1)),
  quan_fig1 %+% subset(quan_df_melt, variable == "i_p_mean" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste( hat(bar(I))[P] )), 
         title = "<b>F.</b> Mature plant induction" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.25), breaks = seq(0,0.2,by = 0.1)),
  quan_fig1 %+% subset(quan_df_melt, variable == "H_SP" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste(hat(H), (hat(S) + hat(P))^-1)), 
         title = "<b>G.</b> Herbivore density" )+
    scale_y_continuous(expand = c(0,0), limits = c(0,6), breaks = seq(0,6,by=2)),
  quan_fig1 %+% subset(quan_df_melt, variable == "P.prop" & theta == 2 & ep_s == 30) +
    labs(y = expression(paste(hat(P), (hat(S) + hat(P))^-1)), 
         title = "<b>H.</b> Prop. mature" )+
    scale_y_continuous(expand = c(0,0), limits = c(0,1),breaks = c(0,0.5,1)),
  nrow = 1)

# Panel row labels.
lab_base <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Baseline", 
           size = 5.5, hjust = 0, family = "sans") +
  theme_void() +
  theme( plot.margin = margin(0,0,0,-750) )

lab_theta <- ggplot() +
  annotate("text", x = 0, y = 0, 
           label = expression(paste("Increased reproduction-defense tradeoff (", theta, " = 4)")), 
           parse = TRUE, size = 5.5, hjust = 0, family = "sans") +
  theme_void() +
  theme( plot.margin = margin(0,0,0,-750) )

Fig3 <- grid.arrange(lab_base, quan_base, 
                     quan_both,
                     lab_theta, quan_theta,
                     nrow = 5, heights = c(0.15,1,1,0.15,1.2))


## Figure 4. ===================================================================
## H/(S+P) ~ tau_s x tau_p for four combinations of m & g.

# Set theme for figures.
theme_uneven <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_markdown(size = 14, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(0,10,3,3),
  legend.position = "none",
  panel.spacing = unit(0.5, "lines")
) }

# Read in simulation output.
quan_df_tau <- read.csv("09_taus-taup-theta1.5-ep30-quan.csv")

# Select only four combinations of m and g to plot.
quan_df_tau <- subset(quan_df_tau,
                      (m == 0.05 | m == 0.8) & (g == 0.15 | g == 0.5))

# Add column of formatted labels for panel headings.
quan_df_tau$type <- with(quan_df_tau, paste(m, g, sep = "-"))
quan_df_tau$type.lab <- factor(quan_df_tau$type,
                        levels = c("0.05-0.15", "0.05-0.5", "0.8-0.15", "0.8-0.5"),
                        labels = c(
                          expression(paste(bold(A.), " ", italic(m), " = 0.05, ", italic(g), " = 0.15")),
                          expression(paste(bold(B.), " ", italic(m), " = 0.05, ", italic(g), " = 0.50")),
                          expression(paste(bold(C.), " ", italic(m), " = 0.8, ", italic(g), " = 0.15")),
                          expression(paste(bold(D.), " ", italic(m), " = 0.8, ", italic(g), " = 0.50"))
                        )
)

# Plot Fig. 4.
Fig4 <- ggplot(data = quan_df_tau, 
                 aes(x = tau_s, y = tau_p, z = H_SP, fill = H_SP)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(color = after_stat(level)), size = 0.4) +
  metR::geom_text_contour(size = 3, skip = 1, stroke = 0.01,
                          label.placer = metR::label_placer_flattest(ref_angle = -60))+
  theme_uneven() + theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, hjust = 0)
  ) +
  labs(x = "Seedling induction<br>responsiveness&nbsp;(&tau;<sub>S</sub>)", 
       y = "Mature plant induction<br>responsiveness&nbsp;(&tau;<sub>P</sub>)", title = "") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_fill_gradient(high = "blue3", low = "gray95", guide = "none")+
  scale_color_gradient(high = "blue4", low = "gray75", guide = "none") +
  facet_wrap(~type.lab, labeller = "label_parsed")
