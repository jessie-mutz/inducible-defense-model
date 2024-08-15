## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Figures of qualitative behavior
## JM, Last updated 2024-01-13

# Load required packages.
require(ggplot2)
require(gridExtra)
require(ggtext)

## Figure 1.====================================================================
## Panels A-C: Behavior ~ tau x epsilon for baseline, increased theta, and 
##  increased r, at low m/low g (m = 0.05 & g = 0.15).
## Panels D-I: Behavior ~ tau x m for baseline, increased theta and increased r, 
##  at low g (g = 0.15 & b = 6) and moderate g (g = 0.5 & b = 2.43).

## (A) Using qualitative simulation output to define regions of model behavior.

# Load output from qualitative simulations for low mature plant mortality (m)/
#  low seedling maturation rate (g) scenario.
load("./inducible-defense-model-v0/07_qualitative-output-lowm-lowg.RData")

# View qualitative behavior (tau x epsilon) from raw simulation output.
ggplot(data = tau_ep_out[["baseline"]]) +
  geom_tile(aes(x = tau_s, y = ep_s, fill = qual))
ggplot(data = tau_ep_out[["theta = 4"]]) +
  geom_tile(aes(x = tau_s, y = ep_s, fill = qual))
ggplot(data = tau_ep_out[["r = 8"]]) +
  geom_tile(aes(x = tau_s, y = ep_s, fill = qual))

# Making it pretty: dataframe for boundaries of model behavior.
tau.range <- unique(tau_ep_out[[1]][,'tau_s'])
epsilon_bounds <- data.frame( cbind( tau = tau.range, base_ext_max = NA, 
                                     theta_pt_min = NA, theta_cyc_min = NA, 
                                     r_ext_max = NA, r_cyc_min = NA ))

# For baseline, increased theta, and increased r parameter sets, find values of 
#  epsilon that defines boundaries of extinction, stability, and cycling at each
#  value of tau.
for ( i in 1:length(tau.range) ) {
  
  # Boundaries of model behavior (tau x epsilon) for baseline parameter set
  temp.df <- subset(tau_ep_out[["baseline"]], tau_s == tau.range[i])
  epsilon_bounds$base_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.25
  
  # Boundaries of model behavior (tau x epsilon) for increased theta
  temp.df <- subset(tau_ep_out[["theta = 4"]], tau_s == tau.range[i])
  epsilon_bounds$theta_pt_min[i] <- min( temp.df[which(temp.df$qual == "pos.pt.eq"), 'ep_s'] ) - 0.25
  if("cycles" %in% temp.df$qual) { 
    epsilon_bounds$theta_cyc_min[i] <- min( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) - 0.25
  }
  
  # Boundaries of model behavior (tau x epsilon) for increased r
  temp.df <- subset(tau_ep_out[["r = 8"]], tau_s == tau.range[i])
  epsilon_bounds$r_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.25
  if("cycles" %in% temp.df$qual) { 
    epsilon_bounds$r_cyc_min[i] <- min( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) - 0.25
  }
  
}

# Use 'loess' to fit smoothed boundary lines for plotting.
epsilon_bounds[sapply(epsilon_bounds, is.infinite)] <- NA

for (i in 2:length(colnames(epsilon_bounds))) {
  col.id <- colnames(epsilon_bounds)[i]
  epsilon_bounds[is.na(epsilon_bounds[,col.id]) == F, paste0(col.id, "_fit")] <- loess(formula = epsilon_bounds[,col.id] ~ tau, data = epsilon_bounds)$fitted
}

# View qualitative behavior (tau x m) from raw simulation output.
ggplot(data = tau_m_out[["baseline"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))
ggplot(data = tau_m_out[["theta = 4"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))
ggplot(data = tau_m_out[["r = 8"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))

# Making it pretty: dataframe for boundaries of model behavior.
tau.range <- unique(tau_m_out[[1]]$tau_s)
m_bounds_lowg <- data.frame(cbind( tau = tau.range, base_ext_min = NA,
                                   theta_ext_min = NA, theta_cyc_max = NA, 
                                   r_ext_min = NA, r_cyc_max = NA ))

# For baseline, increased theta, and increased r parameter sets, find value of 
#  m that defines boundaries of model behavior at each value of tau.
for ( i in 1:length(tau.range) ) {
  
  # Boundaries of model behavior (tau x m) for baseline parameter set
  temp.df <- subset(tau_m_out[["baseline"]], tau_s == tau.range[i])
  m_bounds_lowg$base_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  
  # Boundaries of model behavior (tau x m) for increased theta
  temp.df <- subset(tau_m_out[["theta = 4"]], tau_s == tau.range[i])
  if(tau.range[i] <= 0.42) {
    m_bounds_lowg$theta_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  }
  if("cycles" %in% temp.df$qual) { 
    m_bounds_lowg$theta_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'm'] ) + 0.01
  }
  
  # Boundaries of model behavior (tau x m) for increased r
  temp.df <- subset(tau_m_out[["r = 8"]], tau_s == tau.range[i])
  if("extinct" %in% temp.df$qual) {
    m_bounds_lowg$r_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  }
  if("cycles" %in% temp.df$qual) { 
    m_bounds_lowg$r_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'm'] ) + 0.01
  }
  
}

# Use 'loess' to fit smoothed boundary lines for plotting.
m_bounds_lowg[sapply(m_bounds_lowg, is.infinite)] <- NA

for (i in 2:length(colnames(m_bounds_lowg))) {
  col.id <- colnames(m_bounds_lowg)[i]
  m_bounds_lowg[is.na(m_bounds_lowg[,col.id]) == F, paste0(col.id, "_fit")] <- loess(formula = m_bounds_lowg[,col.id] ~ tau, data = m_bounds_lowg)$fitted
}

# Load output from qualitative simulations for low mature plant mortality (m)/
#  moderate seedling maturation rate (g) scenario.
load("./inducible-defense-model-v0/09_qualitative-output-lowm-modg.RData")

# View qualitative behavior (tau x m) from raw simulation output.
ggplot(data = tau_m_out[["baseline"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))
ggplot(data = tau_m_out[["theta = 4"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))
ggplot(data = tau_m_out[["r = 8"]]) +
  geom_tile(aes(x = tau_s, y = m, fill = qual))

# Making it pretty: dataframe for boundaries of model behavior.
tau.range <- unique(tau_m_out[[1]]$tau_s)
m_bounds_modg <- data.frame(cbind( tau = tau.range, base_ext_min = NA,
                                   theta_ext_min = NA, theta_cyc_min = NA, 
                                   r_ext_min = NA))

# For baseline, increased theta, and increased r parameter sets, find value of 
#  m that defines boundaries of model behavior at each value of tau.
for ( i in 1:length(tau.range) ) {
  
  # Boundaries of model behavior (tau x m) for baseline parameter set
  temp.df <- subset(tau_m_out[["baseline"]], tau_s == tau.range[i])
  m_bounds_modg$base_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  
  # Boundaries of model behavior (tau x m) for increased theta
  temp.df <- subset(tau_m_out[["theta = 4"]], tau_s == tau.range[i])
  m_bounds_modg$theta_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  if("cycles" %in% temp.df$qual) { 
    m_bounds_modg$theta_cyc_min[i] <- min( temp.df[which(temp.df$qual == "cycles"), 'm'] ) - 0.01
  }
  
  # Boundaries of model behavior (tau x m) for increased r
  temp.df <- subset(tau_m_out[["r = 8"]], tau_s == tau.range[i])
  if("extinct" %in% temp.df$qual) {
    m_bounds_modg$r_ext_min[i] <- min( temp.df[which(temp.df$qual == "extinct"), 'm'] ) - 0.01
  }
  
}

# Use 'loess' to fit smoothed boundary lines for plotting.
m_bounds_modg[sapply(m_bounds_modg, is.infinite)] <- NA

for (i in 2:length(colnames(m_bounds_modg))) {
  col.id <- colnames(m_bounds_modg)[i]
  m_bounds_modg[is.na(m_bounds_modg[,col.id]) == F, paste0(col.id, "_fit")] <- loess(formula = m_bounds_modg[,col.id] ~ tau, data = m_bounds_modg)$fitted
}

## (B) Plotting.

# Set theme for figures.
theme_qual <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 14, color = "black"),
  axis.title = element_text(size = 16, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(2,10,2,2)
) }

# Setup for Panels A-C (tau x epsilon).
ep_figs <- ggplot(data = epsilon_bounds) +
  geom_ribbon( aes( x = tau, ymin = 5, ymax = 30 ), fill = "white") +
  theme_qual() +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(5,30), breaks = seq(5,30,by = 5)) +
  labs(x = expression(tau), y = expression(epsilon))

# Setup for Panels D-I (tau x m).
m_figs <- ggplot() +
  geom_ribbon( aes( x = tau, ymin = 0, ymax = 1 ), fill = "white") +
  theme_qual() +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
  labs(x = expression(tau), y = "m")

# Setup for column headings.
lab_figs <- ggplot() + 
  theme_void() +
  scale_y_continuous(limits = c(20,30)) +
  scale_x_continuous(limits = c(0,1)) +
  theme( plot.margin = margin(0,0,0,0) )

# Plot everything.
Fig1 <- grid.arrange(

  # Baseline column heading
  lab_figs + 
    annotate("text", x = 0.1, y = 25, label = "Baseline", size = 5, 
             hjust = 0, family = "serif"),
  
  # Increased theta column heading
  lab_figs + 
    annotate("text", x = 0.1, y = 25, label = "Increased reproduction-", 
             size = 5, hjust = 0, family = "serif") +
    annotate("text", x = 0.1, y = 22, 
             label = expression(paste("defense tradeoff (", theta, " = 4)")), 
             parse = TRUE, size = 5, hjust = 0, family = "serif"), 
  
  # Increased r column heading
  lab_figs +
    annotate("text", x = 0.1, y = 25, 
             label = "Increased herbivore intrinsic", 
             size = 5, hjust = 0, family = "serif") +
    annotate("text", x = 0.1, y = 22, 
             label = expression(paste("growth rate (", italic(r), " = 8)")), 
             parse = TRUE, size = 5, hjust = 0, family = "serif"),   

# Plot and label regions for tau x epsilon using smoothed boundaries in 
# 'epsilon_bounds' dataframe
# Panel A: tau x epsilon, baseline
ep_figs +
  geom_ribbon( aes( x = tau, ymin = 5, ymax = base_ext_max_fit ), fill = "gray65") +
  annotate("text", x = 0.75, y = 6.5, label = "extinct", size = 5, fontface = "italic") +
  annotate("text", x = 0.5, y = 20, label = "stable", size = 5, fontface = "italic") +
  labs(title = "<b>A.</b> <i>g</i> = 0.15, <i>b</i> = 6, <i>m</i> = 0.05"),

#Panel B: tau x epsilon, increased theta
ep_figs +
  geom_ribbon( aes( x = tau, ymin = 5, ymax = theta_pt_min_fit ), fill = "gray65") +
  geom_ribbon( aes( x = tau, ymin = theta_cyc_min_fit, ymax = theta_pt_min_fit), fill = "gray85") +
  annotate("text", x = 0.5, y = 14, label = "extinct", size = 5, fontface = "italic") +
  annotate("text", x = 0.75, y = 26.5, label = "cycles", size = 5, fontface = "italic") +
  annotate("text", x = 0.025, y = 28, label = "stable", size = 5, hjust = 0, fontface = "italic") +
  labs(title = "<b>B.</b> <i>g</i> = 0.15, <i>b</i> = 6, <i>m</i> = 0.05"),

#Panel C: tau x epsilon, increased r
ep_figs +
  geom_ribbon( aes( x = tau, ymin = 5, ymax = r_ext_max_fit ), fill = "gray65") +
  geom_ribbon( aes( x = tau, ymin = r_cyc_min_fit, ymax = 30 ), fill = "gray85") +
  annotate("text", x = 0.75, y = 6.5, label = "extinct", size = 5, fontface = "italic") +
  annotate("text", x = 0.5, y = 16, label = "stable", size = 5, fontface = "italic") +
  annotate("text", x = 0.8, y = 26.5, label = "cycles", size = 5, fontface = "italic") +
  labs(title = "<b>C.</b> <i>g</i> = 0.15, <i>b</i> = 6, <i>m</i> = 0.05"),

# Plot and label regions for tau x m (low maturation scenario) using smoothed 
# boundaries in 'm_bounds_lowg' dataframe
# Panel D: tau x m, baseline
m_figs %+% m_bounds_lowg + 
  geom_ribbon( aes( x = tau, ymin = base_ext_min_fit, ymax = 1 ), fill = "gray65") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.025, y = 0.85, label = "ext.", hjust = 0, size = 5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.65, label = "stable", size = 5, fontface = "italic") +
  labs(title = "<b>D.</b> <i>g</i> = 0.15, <i>b</i> = 6"),

# Panel E: tau x m, increased theta
m_figs %+% m_bounds_lowg + 
  geom_ribbon( aes( x = tau, ymin = 0, ymax = 1 ), fill = "gray65") +
  geom_ribbon( aes( x = tau, ymin = 0, ymax = theta_ext_min_fit ), fill = "white") +
  geom_ribbon( aes( x = tau, ymin = 0, ymax = theta_cyc_max_fit ), fill = "gray85") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.5, y = 0.7, label = "extinct", size = 5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.15, label = "cycles", size = 5, fontface = "italic") +
  annotate("text", x = 0.078, y = 0.11, label = "stable", size = 5, hjust = 0, fontface = "italic") +
  labs(title = "<b>E.</b> <i>g</i> = 0.15, <i>b</i> = 6"),

# Panel F: tau x m, increased r
m_figs %+% m_bounds_lowg + 
  geom_ribbon( aes( x = tau, ymin = r_ext_min_fit, ymax = 1 ), fill = "gray65") +
  geom_ribbon( aes( x = tau, ymin = 0, ymax = r_cyc_max_fit ), fill = "gray85") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.75, y = 0.3, label = "cycles", size = 5, fontface = "italic") +
  annotate("text", x = 0.025, y = 0.85, label = "ext.", hjust = 0, size = 5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.65, label = "stable", size = 5, fontface = "italic") +
  annotate("segment", x = 0.75, xend = 0.75, y = 0.2, yend = 0.1,
           arrow = arrow(length = unit(0.3, "cm"))) +
  labs(title = "<b>F.</b> <i>g</i> = 0.15, <i>b</i> = 6"),

# Plot and label regions for tau x m (moderate maturation scenario) using smoothed 
# boundaries in 'm_bounds_modg' dataframe
# Panel G: tau x m, baseline
m_figs %+% m_bounds_modg + 
  geom_ribbon( aes( x = tau, ymin = base_ext_min_fit, ymax = 1 ), fill = "gray65") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.05, y = 0.85, label = "ext.", hjust = 0, size = 5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.65, label = "stable", size = 5, fontface = "italic") +
  labs(title = "<b>G.</b> <i>g</i> = 0.50, <i>b</i> = 2.43"),

# Panel H: tau x m, increased theta
m_figs %+% m_bounds_modg + 
  geom_ribbon( aes( x = tau, ymin = theta_ext_min_fit, ymax = 1 ), fill = "gray65") +
  geom_ribbon( aes( x = tau, ymin = theta_cyc_min_fit, ymax = theta_ext_min_fit ), fill = "gray85") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.5, y = 0.7, label = "extinct", size = 5, fontface = "italic") +
  annotate("text", x = 0.75, y = 0.35, label = "cycles", size = 5, fontface = "italic") +
  annotate("text", x = 0.25, y = 0.3, label = "stable", size = 5, fontface = "italic") +
  annotate("segment", x = 0.25, xend = 0.25, y = 0.23, yend = 0.12,
           arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("segment", x = 0.75, xend = 0.75, y = 0.28, yend = 0.12,
           arrow = arrow(length = unit(0.3, "cm"))) +
  labs(title = "<b>H.</b> <i>g</i> = 0.50, <i>b</i> = 2.43"),

# Panel I: tau x m, increased r
m_figs %+% m_bounds_modg + 
  geom_ribbon( aes( x = tau, ymin = r_ext_min_fit, ymax = 1 ), fill = "gray65") +
  geom_hline( yintercept = 0.05, linetype = "dashed") +
  annotate("text", x = 0.05, y = 0.85, label = "ext.", hjust = 0, size = 5, fontface = "italic") +
  annotate("text", x = 0.7, y = 0.65, label = "stable", size = 5, fontface = "italic") +
  labs(title = "<b>I.</b> <i>g</i> = 0.50, <i>b</i> = 2.43"),

nrow = 4, heights = c(0.35,1,1,1))


## Figure S1.===================================================================
## Min. epsilon for persistence (for baseline) and area of cycling (for increased 
##  theta) ~ tau for four combinations of m & g.

## (A) Using qualitative simulation output to define boundaries of model behavior.

# Dataframe for boundaries of model behavior.
epsilon_bounds2 <- data.frame( cbind( tau = tau.range, 
                                          base_ext_min = NA,
                                          theta_ext_max = NA,
                                          theta_cyc_max = NA ))

# Make a list containing a dataframe for each life history scenario.
epsilon_bounds2 <- list(epsilon_bounds2, epsilon_bounds2, epsilon_bounds2, epsilon_bounds2)
names(epsilon_bounds2) <- c("lowm_lowg", "highm_lowg", "lowm_modg", "highm_modg")

for ( i in 1:length(tau.range) ) {

  # Boundaries of model behavior for low m & low g
  tau_ep_out <- readRDS("2023-07-28_USDA-cotton_low-m-low-g-qual-ep.RDS")
  temp.df <- subset(tau_ep_out[["baseline"]], tau_s == tau.range[i])
  epsilon_bounds2[[1]]$base_ext_min[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  
  temp.df <- subset(tau_ep_out[["theta = 4"]], tau_s == tau.range[i])
  epsilon_bounds2[[1]]$theta_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  if("cycles" %in% temp.df$qual) {
    epsilon_bounds2[[1]]$theta_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) + 0.5
  }
  
  # Boundaries of model behavior for high m & low g
  tau_ep_out <- readRDS("2023-07-28_USDA-cotton_high-m-low-g-qual-ep.RDS")
  temp.df <- subset(tau_ep_out[["baseline"]], tau_s == tau.range[i])
  epsilon_bounds2[[2]]$base_ext_min[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  
  temp.df <- subset(tau_ep_out[["theta = 4"]], tau_s == tau.range[i])
  epsilon_bounds2[[2]]$theta_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  if("cycles" %in% temp.df$qual) {
    epsilon_bounds2[[2]]$theta_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) + 0.5
  }
  
  # Boundaries of model behavior for low m & moderate g
  tau_ep_out <- readRDS("2023-07-28_USDA-cotton_low-m-mod-g-qual-ep.RDS")
  temp.df <- subset(tau_ep_out[["baseline"]], tau_s == tau.range[i])
  epsilon_bounds2[[3]]$base_ext_min[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  
  temp.df <- subset(tau_ep_out[["theta = 4"]], tau_s == tau.range[i])
  epsilon_bounds2[[3]]$theta_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  if("cycles" %in% temp.df$qual) {
    epsilon_bounds2[[3]]$theta_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) + 0.5
  }
  
  # Boundaries of model behavior for high m & moderate g
  tau_ep_out <- readRDS("2023-07-28_USDA-cotton_high-m-mod-g-qual-ep.RDS")
  temp.df <- subset(tau_ep_out[["baseline"]], tau_s == tau.range[i])
  epsilon_bounds2[[4]]$base_ext_min[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  
  temp.df <- subset(tau_ep_out[["theta = 4"]], tau_s == tau.range[i])
  epsilon_bounds2[[4]]$theta_ext_max[i] <- max( temp.df[which(temp.df$qual == "extinct"), 'ep_s'] ) + 0.5
  if("cycles" %in% temp.df$qual) {
    epsilon_bounds2[[4]]$theta_cyc_max[i] <- max( temp.df[which(temp.df$qual == "cycles"), 'ep_s'] ) + 0.5
  }
  
  }

# Use 'loess' to fit smoothed boundary lines for plotting.
for (j in 1:4) {
  epsilon_bounds2[[j]][sapply(epsilon_bounds2[[j]], is.infinite)] <- NA
for (i in 2:length(colnames(epsilon_bounds2[[j]]))) {
  col.id <- colnames(epsilon_bounds2[[j]])[i]
  epsilon_bounds2[[j]][is.na(epsilon_bounds2[[j]][,col.id]) == F, paste0(col.id, "_fit")] <- loess(formula = epsilon_bounds2[[j]][,col.id] ~ tau, data = epsilon_bounds2[[j]])$fitted
}
}

epsilon_bounds2[[4]]$theta_cyc_max_fit <- NA

# Combine list element into single dataframe.
ep_test <- rbind(cbind(epsilon_bounds2[[1]], m = "low", g = "low"),
                 cbind(epsilon_bounds2[[2]], m = "high", g = "low"),
                 cbind(epsilon_bounds2[[3]], m = "low", g = "mod"),
                 cbind(epsilon_bounds2[[4]], m = "high", g = "mod")
)

## (B) Plotting.

# Plot min. value of epsilon needed for persistence for each life history 
#  scenario(Fig. S1A).
ep1 <- ggplot( data = ep_test ) +
  geom_line( aes( x = tau, y = base_ext_min_fit, color = g, linetype = m)) +
  annotate("text", x = 0.65, y = 35, 
           label = expression(paste(italic(m), " = 0.05")), hjust = 0) +
  annotate("text", x = 0.65, y = 31, 
           label = expression(paste(italic(m), " = 0.80")), hjust = 0) +
  annotate("text", x = 0.65, y = 27, 
           label = expression(paste(italic(g), " = 0.15")), hjust = 0) +
  annotate("text", x = 0.65, y = 23, 
           label = expression(paste(italic(g), " = 0.50")), hjust = 0) +
  annotate("text", x = 0.08, y = 3.5, label = "extinct", size = 5, hjust = 0, fontface = "italic") +
  annotate("text", x = 0.08, y = 18, label = "stable", size = 5, hjust = 0, fontface = "italic") +
  annotate("segment", x = 0.5, xend = 0.6, y = 35, yend = 35) +
  annotate("segment", x = 0.5, xend = 0.6, y = 31, yend = 31, linetype = "dashed") +
  annotate("point", x = 0.55, y = 27, shape = 22, fill = "gray80", color = "white", size = 5.5) +
  annotate("point", x = 0.55, y = 23, shape = 22, fill = "black", color = "white", size = 5.5) +
  theme_qual() + theme(
    plot.margin = margin(5,12,5,2)
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = seq(0, 40, by = 10)) +
  labs(x = expression(tau), y = expression(paste(epsilon)),
       title = "<b>A.</b> Baseline") +
  scale_color_manual(values = c("gray", "black"), breaks = c("low", "mod"), guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed"), breaks = c("low", "high"), guide = "none")

# Plot regions of cycling for each life history scenario (Fig. S1B).
ep2 <- ggplot( data = ep_test ) +
  geom_line( aes( x = tau, y = theta_ext_max_fit, color = g, linetype = m)) +
  geom_line( aes( x = tau, y = theta_cyc_max_fit, color = g, linetype = m)) +

  geom_ribbon( aes( x = tau, ymin = theta_ext_max_fit, ymax = theta_cyc_max_fit,
                    fill = g, alpha = m)) +
  annotate("segment", xend = 0.7, x = 0.6, yend = 25.5, y = 16,
           arrow = arrow(length = unit(0.5,"lines")), size = 0.6 ) +
  annotate("segment", xend = 0.5, x = 0.6, yend = 24, y = 16,
           arrow = arrow(length = unit(0.5,"lines")), size = 0.6 ) +
  annotate("segment", xend = 0.9, x = 0.6, yend = 28, y = 16,
           arrow = arrow(length = unit(0.5,"lines")), size = 0.6 ) +
  annotate("text", x = 0.6, y = 14, label = "cycles", size = 5, fontface = "italic") +
  annotate("text", x = 0.08, y = 5, label = "extinct", size = 5, hjust = 0, fontface = "italic") +
  annotate("text", x = 0.08, y = 35, label = "stable", size = 5, hjust = 0, fontface = "italic") +
  theme_qual() + theme(
    plot.margin = margin(5,12,5,2)
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = c(0,0.5,1)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = seq(0,40,by = 10)) +
  labs(x = expression(tau), y = expression(epsilon),
       title = "<b>B.</b> &theta; = 4") +
  scale_color_manual(values = c("gray", "black"), breaks = c("low", "mod"), guide = "none") +
  scale_linetype_manual(values = c("solid", "dashed"), breaks = c("low", "high"), guide = "none") +
  scale_fill_manual(values = c("gray", "black"), breaks = c("low", "mod"), guide = "none") +
  scale_alpha_manual(values = c(0.3, 0.1), breaks = c("low", "high"), guide = "none")

FigS1 <- grid.arrange(ep1, ep2, ncol = 1)
