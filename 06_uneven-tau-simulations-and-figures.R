## Life history modulates effects of inducible defenses on consumer-resource dynamics.
## Supplemental code: Simulations and figures of uneven tau
## JM, Last updated 2024-01-04

# Source model functions for model simulations.
source("./inducible-defense-model/01_model-functions.R")

# Load required packages.
require(reshape2)
require(doBy)
require(ggplot2)
require(metR)
require(ggtext)

# Set theme for figures.
theme_uneven <- function() { theme_bw() + theme(
  panel.grid = element_blank(),
  axis.text = element_text(size = 12, color = "black"),
  axis.title = element_text(size = 14, family = "serif"),
  plot.title = element_markdown(size = 12.5),
  plot.margin = margin(0,10,3,3),
  legend.position = "none",
  panel.spacing = unit(0.5, "lines")
) }

## Simulations for the below figures.===========================================

# Initial parameter list.
params <- c(
  b = 6,
  theta = 1.5,
  g = 0.15,
  x = 1,
  tau_s = 1,
  tau_p = 1,
  c = 0.8,
  a = 0.1,
  ep_s = 27,
  ep_p = 27,
  w = 1,
  r = 2,
  mu = 0.15,
  m = 0.05
)

times <- seq(0,2500, by = 2)

# Set list of parameters to vary.
param.vals <- list(
  'ep_s' = 30,
  'ep_s' = 30,
  'tau_s' = seq(0, 1, by = 0.05),
  'tau_p' = seq(0, 1, by = 0.05)
)

# List with combinations of b, m, and g to use.
list1 <- list(
  c(6, 0.05, 0.15),
  c(11, 0.8, 0.15),
  c(2.43, 0.05, 0.5),
  c(4.455, 0.8, 0.5)
)

# For each element of list1, reassign b, m, and g, find equilibrium as a 
#  function of tau_s and tau_p, and save output to an object.

for (i in 1:length(list1)) {
  
  params[c('b', 'm', 'g')] <- list1[[i]]
  
  out1 <- funcy('tau_s', 'tau_p', c('ep_s', 'ep_p') )
  out2 <- reshape2::dcast(out1, tau_s + tau_p + ep_s + timestep ~ variable, margin = "value")
 
  # Add columns for prop. of seedling biomass (S.prop), mean induction of
  #  seedlings (i_s_mean) and mature plants (i_p_mean), and herbivore density per
  #  unit plant biomass (H_SP).
  out2[c('S.prop', 'i_s_mean', 'i_p_mean', 'H_SP')] <- c( 
    with(out2, S / (S + P)),
    with(out2, I_s / S),
    with(out2, I_p / P),
    with(out2, H / (S + P)))

  out2.sum <- doBy::summaryBy(S + P + i_s_mean + i_p_mean + H_SP + H ~ tau_s + tau_p,
                              data = out2, FUN = mean, args = c(na.rm = T))

  assign(paste0("tstp_m", params['m'], '_g', params['g']), out2.sum)
}

# Combine objects into a single dataframe.
comb <- rbind( cbind(tstp_m0.05_g0.15, m = "low", g = "low"),
               cbind(tstp_m0.8_g0.15, m = "high", g = "low"),
               cbind(tstp_m0.05_g0.5, m = "low", g = "mod"),
               cbind(tstp_m0.8_g0.5, m = "high", g = "mod")
)
comb$type <- paste0(comb$m, "_", comb$g)
comb$type <- factor(comb$type, levels = c("low_low", "low_mod", "high_low", "high_mod"))

# Interpolate to replace one NA value.
comb[which(comb$tau_s == 0 & comb$tau_p == 0.55 & comb$type == "low_mod"),
     c('S.mean', 'P.mean', 'i_s_mean.mean', 'i_p_mean.mean', 'H_SP.mean', 'H.mean')] <-
  c(loess(formula = S.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12],
loess(formula = P.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12],
loess(formula = i_s_mean.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12],
loess(formula = i_p_mean.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12],
loess(formula = H_SP.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12],
loess(formula = H.mean ~ tau_p, data = comb[which(comb$tau_s == 0 & comb$type == "low_mod"),])$fitted[12])


## Figure 5. =======================================================
## H/(S+P) ~ tau_s x tau_p for four combinations of m & g.

# Add column of formatted labels for panel headings.
comb$type.lab <- factor(comb$type,
                        levels = c("low_low", "low_mod", "high_low", "high_mod"),
                        labels = c(
                          expression(paste(bold(A.), " ", italic(m), " = 0.05, ", italic(g), " = 0.15")),
                          expression(paste(bold(B.), " ", italic(m), " = 0.05, ", italic(g), " = 0.50")),
                          expression(paste(bold(C.), " ", italic(m), " = 0.8, ", italic(g), " = 0.15")),
                          expression(paste(bold(D.), " ", italic(m), " = 0.8, ", italic(g), " = 0.50"))
                        )
)

# Plot Fig. 5.
facet1 <- ggplot(data = comb, aes(x = tau_s, y = tau_p, z = H_SP.mean, fill = H_SP.mean)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(color = ..level..), size = 0.4) +
  metR::geom_text_contour(size = 3, skip = 1, stroke = 0.01,
                          label.placer = label_placer_flattest(ref_angle = -60))+
  theme_uneven() + theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, hjust = 0)
  ) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_fill_gradient(high = "blue3", low = "gray95", guide = "none")+
  scale_color_gradient(high = "blue4", low = "gray75", guide = "none") +
  facet_wrap(~type.lab, labeller = "label_parsed")


## Figure S8. =================================================================
## State variables at equilibrium ~ tau_s x tau_p for four combinations of m & g.

# New facet label names.
type.labs <-  c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", 
                "I.", "J.", "K.", "L.", "M.", "N.", "O.", "P.",
                "Q.", "R.", "S.", "T.", "U.", "V.", "W.", "X.")
names(type.labs) <- rep(c("low_low", "low_mod", "high_low", "high_mod"), 6)

# (uses 'facet1' plotting environment defined in Fig. 5 above)
FigS8 <- grid.arrange(

# Seedling biomass as a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = S.mean, fill = S.mean) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_gradient(high = "blue3", low = "gray95", guide = "none", limits = c(0,5)) +
  scale_color_gradient(high = "blue4", low = "gray75", guide = "none", limits = c(0,5)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[1:4]), nrow = 1),
  
# Mature plant biomass as a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = P.mean, fill = P.mean) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_gradient(high = "blue3", low = "white", guide = "none", limits = c(0,1)) +
  scale_color_gradient(high = "blue4", low = "gray90", guide = "none", limits = c(0,1)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[5:8]), nrow = 1),

# Mean induction of seedlings a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = i_s_mean.mean, fill = i_s_mean.mean) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_gradient(high = "blue3", low = "white", guide = "none", limits = c(-1e-4,0.22)) +
  scale_color_gradient(high = "blue4", low = "gray90", guide = "none", limits = c(-1e-4,0.22)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[9:12]), nrow = 1),

# Mean induction of mature plants as a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = i_p_mean.mean, fill = i_p_mean.mean) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_gradient(high = "blue3", low = "white", guide = "none", limits = c(-1e-4,0.4)) +
  scale_color_gradient(high = "blue4", low = "gray90", guide = "none", limits = c(-1e-4,0.4)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[13:16]), nrow = 1),

# Herbivore population size as a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = H.mean, fill = H.mean) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_gradient(high = "blue3", low = "white", guide = "none", limits = c(0,14)) +
  scale_color_gradient(high = "blue4", low = "gray90", guide = "none", limits = c(0,14)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[17:20]), nrow = 1),

# Prop. of mature plant biomass as a function of tau_s and tau_p for each plant life history scenario.
facet1 + aes(x = tau_s, y = tau_p, z = P.mean/(S.mean+P.mean), fill = P.mean/(S.mean+P.mean)) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  theme(plot.margin = margin(2,2,2,2)) +
  scale_fill_gradient(high = "blue3", low = "white", guide = "none", limits = c(0,1)) +
  scale_color_gradient(high = "blue4", low = "gray90", guide = "none", limits = c(0,1)) +
  facet_wrap(~type, labeller = labeller(type = type.labs[21:24]), nrow = 1),

ncol = 1, heights = c(rep(0.85, 5), 1.05))


## Figure S9. =================================================================
## Relative consumption rates at equilibrium ~ tau_s x tau_p for four 
## combinations of m & g.

FigS9 <- ggplot(data = comb, aes(x = tau_s, y = tau_p, 
                        z = ((1 + 30*i_s_mean.mean) + H_SP.mean)/((1 + 30*i_p_mean.mean) + H_SP.mean), 
                        fill = ((1 + 30*i_s_mean.mean) + H_SP.mean)/((1 + 30*i_p_mean.mean) + H_SP.mean) )) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(color = ..level..), size = 0.4) +
  metR::geom_text_contour(size = 3.3, skip = 1, stroke = 0.05,
                          label.placer = label_placer_flattest(ref_angle = -60))+
  theme_uneven() + theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10.5, hjust = 0)
  ) +
  labs(x = expression(tau[S]), y = expression(tau[P]), title = "") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_x_continuous(expand = c(0,0), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_fill_gradient2(high = "blue", low = "gray", mid = "white", midpoint = 1, guide = "none") +
  scale_color_gradient2(high = "blue", low = "gray40", mid = "gray80", midpoint = 1, guide = "none") +
  facet_wrap(~type.lab, labeller = "label_parsed")

