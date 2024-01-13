## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Simulations and figures of quantitative equilibrium dynamics
## JM, Last updated 2024-01-04

# Source model functions for model simulations.
source("./inducible-defense-model/01_model-functions.R")

# Load required packages.
require(ggplot2)
require(gridExtra)
require(reshape2)

# Set theme for figures.
theme_quan <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 14, color = "black"),
  axis.title = element_text(size = 16, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(3,15,3,3),
  panel.border = element_blank(), 
  axis.line = element_line(),
  legend.position = "none"
) }


## Simulations for the below figures.===========================================

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
  ep_s = 30,
  ep_p = 30,
  w = 1,
  r = 2,
  mu = 0.15,
  m = 0.05
)

# Values of theta, tau, and epsilon to use in simulations.
param.vals <- list(
  'theta' = c(1.5, 4),
  'tau_s' = seq(0, 1, by = 0.05),
  'tau_p' = seq(0, 1, by = 0.05),
  'ep_s' = c(15, 30, 45),
  'ep_p' = c(15, 30, 45)
)

times <- seq(0, 2000, by = 1)

# List with combinations of b, m, and g to use.
list1 <- list(
  c(6, 0.05, 0.15),
  c(11, 0.8, 0.15),
  c(2.43, 0.05, 0.5),
  c(4.455, 0.8, 0.5)
)

# For each element of list1, reassign b, m, and g, find equilibrium as a 
#  function of theta, tau, and epsilon, and save output to an object.

for (i in 1:length(list1)) {
  
  params[c('b', 'm', 'g')] <- list1[[i]]
  
  out1 <- funcy('theta', c('tau_s', 'tau_p'), c('ep_s', 'ep_p') )
  out2 <- reshape2::dcast(out1, theta + tau_s + ep_s + timestep ~ variable, margin = "value")
  
  # Add columns for prop. of mature plant biomass (P.prop), mean induction of
  #  seedlings (i_s_mean) and mature plants (i_p_mean), herbivore density per
  #  unit plant biomass (H_SP), and total no. of seeds produced (seeds).
  out2[c('P.prop', 'i_s_mean', 'i_p_mean', 'H_SP', 'seeds')] <- c( 
    with(out2, P / (S + P)),
    with(out2, I_s / S),
    with(out2, I_p / P),
    with(out2, H / (S + P)),
    with(out2, params['b']*(P - I_p*theta)))
  out2$seeds <- sapply( out2$seeds, function(x) { max(0, x) } )

  out3 <- reshape2::melt(out2, measure.vars = c("S", "P", "P.prop", "i_s_mean", 
                                                "i_p_mean",  "H", "H_SP", "seeds"))
  out3$m <- params['m']
  out3$g <- params['g']
  
  assign(paste0('quan_m', params['m'], '_g', params['g']), out3)
  
}

# Combine objects into a single dataframe.
out3 <- rbind(quan_m0.05_g0.15, quan_m0.8_g0.15, quan_m0.05_g0.5, quan_m0.8_g0.5)

out3[which(out3$g == 0.5 & out3$m == 0.8 & out3$tau_s > 0.5 & out3$theta == 4),'value'] <- NA


## Figure 4.====================================================================

# Set up environment for plotting equilibrium as a function of tau, with lines
#  for each combination of m and g values.
quan_fig1 <- ggplot( ) +
  geom_line(aes( x = tau_s, y = value, color = as.factor(g), linetype = as.factor(m) ))+
  theme_quan() + 
  labs(x = expression(tau), y = expression(paste( hat(S) )), 
       title = "<b>A.</b> Seedling biomass")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,5)) +
  scale_color_manual(values = c("gray", "black")) +
  scale_linetype_manual(values = c("solid", "dashed")) 

# Plotting for baseline theta.
quan_base <- grid.arrange( 
  quan_fig1 %+% subset(out3, variable == "S" & theta == 1.5 & ep_s == 30) +
    annotate("segment", x = 0.5, xend = 0.6, y = 4, yend = 4) +
    annotate("segment", x = 0.5, xend = 0.6, y = 3.5, yend = 3.5, linetype = "dashed") +
    annotate("point", x = 0.55, y = 3, shape = 22, fill = "gray80", color = "white", size = 5.5) +
    annotate("point", x = 0.55, y = 2.5, shape = 22, fill = "black", color = "white", size = 5.5) +  
    annotate("text", x = 0.63, y = 4, label = expression(paste(italic(m), " = 0.05")), hjust = 0) +
    annotate("text", x = 0.63, y = 3.5, label = expression(paste(italic(m), " = 0.80")), hjust = 0) +
    annotate("text", x = 0.63, y = 3, label = expression(paste(italic(g), " = 0.15")), hjust = 0) +
    annotate("text", x = 0.63, y = 2.5, label = expression(paste(italic(g), " = 0.50")), hjust = 0),
  quan_fig1 %+% subset(out3, variable == "P" & theta == 1.5 & ep_s == 30) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )), 
         title = "<b>B.</b> Mature plant biomass"),
  quan_fig1 %+% subset(out3, variable == "H" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), 
         title = "<b>C.</b> Herbivore pop. size" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  quan_fig1 %+% subset(out3, variable == "seeds" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste("b", hat(P), "(1 - ", theta, hat(bar(I))[P],")")), 
         title = "<b>D.</b> Seed production" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,7)),
  nrow = 1)

# Plotting for increased theta.
quan_theta <- grid.arrange( 
  quan_fig1 %+% subset(out3, variable == "S" & theta == 4 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste( hat(S) )), 
         title = "<b>I.</b> Seedling biomass"), 
  quan_fig1 %+% subset(out3, variable == "P" & theta == 4 & ep_s == 30) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )), 
         title = "<b>J.</b> Mature plant biomass"),
  quan_fig1 %+% subset(out3, variable == "H" & theta == 4 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), 
         title = "<b>K.</b> Herbivore pop. size" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  quan_fig1 %+% subset(out3, variable == "seeds" & theta == 4 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste("b", hat(P), "(1 - ", theta, hat(bar(I))[P],")")), 
         title = "<b>L.</b> Seed production" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,7)),
  nrow = 1)

# Plotting for variables with no effect of theta.
quan_both <- grid.arrange(
  quan_fig1 %+% subset(out3, variable == "i_s_mean" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[S] )), 
         title = "<b>E.</b> Seedling induction" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.25), breaks = seq(0,0.2,by = 0.1)),
  quan_fig1 %+% subset(out3, variable == "i_p_mean" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[P] )), 
         title = "<b>F.</b> Mature plant induction" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.25), breaks = seq(0,0.2,by = 0.1)),
  quan_fig1 %+% subset(out3, variable == "H_SP" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste(hat(H), (hat(S) + hat(P))^-1)), 
         title = "<b>G.</b> Herbivore density" )+
    scale_y_continuous(expand = c(0,0), limits = c(0,6), breaks = seq(0,6,by=2)),
  quan_fig1 %+% subset(out3, variable == "P.prop" & theta == 1.5 & ep_s == 30) +
    labs(x = expression(tau), y = expression(paste(hat(P), (hat(S) + hat(P))^-1)), 
         title = "<b>H.</b> Prop. mature" )+
    scale_y_continuous(expand = c(0,0), limits = c(0,1),breaks = c(0,0.5,1)),
  nrow = 1)

# Panel row labels.
lab_base <- ggplot() +
  annotate("text", x = 0, y = 0, label = "Baseline", 
           size = 5, hjust = 0, family = "serif") +
  theme_void() +
  theme( plot.margin = margin(0,0,0,-750) )

lab_theta <- ggplot() +
  annotate("text", x = 0, y = 0, 
           label = expression(paste("Increased reproduction-defense tradeoff (", theta, " = 4)")), 
           parse = TRUE, size = 5, hjust = 0, family = "serif") +
  theme_void() +
  theme( plot.margin = margin(0,0,0,-750) )

Fig4 <- grid.arrange(lab_base, quan_base, 
                     quan_both,
                     lab_theta, quan_theta,
                     nrow = 5, heights = c(0.15,1,1,0.15,1))


## Figure S6.===================================================================
## State variables ~ tau & epsilon at four combinations of m and g.

# Set up environment for plotting equilibrium as a function of tau, with lines
#  for each of three values of epsilon.
quan_fig2 <- ggplot( ) +
  geom_line(aes( x = tau_s, y = value, linetype = as.factor(ep_s) ))+
  theme_quan() +
  labs(x = expression(tau), y = expression(paste( hat(S) )), 
       title = "<b>A.</b> Seedling biomass")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.5)) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))

FigS6 <- grid.arrange( 
  
  # For m = 0.05, g = 0.15
  quan_fig2 %+% subset(out3, variable == "S" & theta == 1.5 & m == 0.05 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(S) )), 
         title = "<b>A.</b> <i>m</i> = 0.05, <i>g</i> = 0.15") +
    annotate("segment", x = 0.5, xend = 0.6, y = 4, yend = 4, linetype = "dotted") +
    annotate("segment", x = 0.5, xend = 0.6, y = 3.2, yend = 3.2, linetype = "dashed") +
    annotate("segment", x = 0.5, xend = 0.6, y = 2.4, yend = 2.4, linetype = "solid") +
    annotate("text", x = 0.63, y = 4, label = expression(paste(epsilon, " = 15")), hjust = 0) +
    annotate("text", x = 0.63, y = 3.2, label = expression(paste(epsilon, " = 30")), hjust = 0) +
    annotate("text", x = 0.63, y = 2.4, label = expression(paste(epsilon, " = 45")), hjust = 0) +
    scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0,5)),
  quan_fig2 %+% subset(out3, variable == "P" & theta == 1.5 & m == 0.05 & g == 0.15) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )), title = "<b>B.</b>"),
  quan_fig2 %+% subset(out3, variable == "i_s_mean" & theta == 1.5 & m == 0.05 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[S] )), title = "<b>C.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "i_p_mean" & theta == 1.5 & m == 0.05 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[P] )), title = "<b>D.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "H" & theta == 1.5 & m == 0.05 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), title = "<b>E.</b>" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  
  # For m = 0.8, g = 0.15
  quan_fig2 %+% subset(out3, variable == "S" & theta == 1.5 & m == 0.8 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(S) )), 
         title = "<b>F.</b> <i>m</i> = 0.80, <i>g</i> = 0.15") +
    scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0,5)),
  quan_fig2 %+% subset(out3, variable == "P" & theta == 1.5 & m == 0.8 & g == 0.15) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )), title = "<b>G.</b>"),
  quan_fig2 %+% subset(out3, variable == "i_s_mean" & theta == 1.5 & m == 0.8 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[S] )), title = "<b>H.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "i_p_mean" & theta == 1.5 & m == 0.8 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[P] )), title = "<b>I.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "H" & theta == 1.5 & m == 0.8 & g == 0.15) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), title = "<b>J.</b>" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  
  # For m = 0.05, g = 0.5
  quan_fig2 %+% subset(out3, variable == "S" & theta == 1.5 & m == 0.05 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(S) )), 
         title = "<b>K.</b> <i>m</i> = 0.05, <i>g</i> = 0.50") +
    scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0,5)),
  quan_fig2 %+% subset(out3, variable == "P" & theta == 1.5 & m == 0.05 & g == 0.5) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )), title = "<b>L.</b>"),
  quan_fig2 %+% subset(out3, variable == "i_s_mean" & theta == 1.5 & m == 0.05 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[S] )), title = "<b>M.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "i_p_mean" & theta == 1.5 & m == 0.05 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[P] )), title = "<b>N.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "H" & theta == 1.5 & m == 0.05 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), title = "<b>O.</b>" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  
  # For m = 0.8, g = 0.5
  quan_fig2 %+% subset(out3, variable == "S" & theta == 1.5 & m == 0.8 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(S) )), 
         title = "<b>P.</b> <i>m</i> = 0.80, <i>g</i> = 0.50") +
    scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0,5)),
  quan_fig2 %+% subset(out3, variable == "P" & theta == 1.5 & m == 0.8 & g == 0.5) +
    scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1))+
    labs(x = expression(tau), y = expression(paste( hat(P) )),  title = "<b>Q.</b>"),
  quan_fig2 %+% subset(out3, variable == "i_s_mean" & theta == 1.5 & m == 0.8 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[S] )), title = "<b>R.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "i_p_mean" & theta == 1.5 & m == 0.8 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(bar(I))[P] )), title = "<b>S.</b>" )+
    scale_y_continuous(expand = c(0,0), limits = c(-1e-4,0.45), breaks = seq(0,0.4,by = 0.2)),
  quan_fig2 %+% subset(out3, variable == "H" & theta == 1.5 & m == 0.8 & g == 0.5) +
    labs(x = expression(tau), y = expression(paste( hat(H) )), title = "<b>T.</b>" ) +
    scale_y_continuous(expand = c(0,0), limits = c(0,14)),
  
  nrow = 4)
