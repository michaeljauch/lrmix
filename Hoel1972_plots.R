#############
# Hoel data #
#############
# IMPORTANT:
# Download the following RData file:
# https://www.dropbox.com/s/z4jcqfaf49tb342/sample_densities-1.RData?dl=0
# Then, load it with load()
# Alternatively, one could run the MCMC algorithm by
# running the script Hoel1972_runMCMC.R 

library(tidyverse)
library(ggthemes)
library(gridExtra)


# X: control group (f)
# Y: germ-free group (g)
# Assuming X <= Y in the LR order  

df_XY = data.frame(x = c(X, Y),
                   group = rep(c("control", "germ-free"), c(length(X), length(Y))))

# plot data
ggplot(df_XY) +
  aes(x = x) +
  geom_histogram(alpha = 0.5) +
  theme_tufte() +
  facet_wrap(~group, nrow = 2)


#########
# theta #
#########
theta = sapply(pred, with, theta)
mean(theta == 1) # prop of thetas equal to 1
df_theta = data.frame(theta = theta)
summary(theta)
# plotting theta
ggplot(df_theta) +
  aes(x = theta) +
  geom_histogram(aes(y=..density..), alpha = 0.5) +
  geom_density(alpha=.2, lwd = 1) +
  theme_tufte() +
  xlab(expression(theta))

#################
# ratio r = f/g #
#################
f = sapply(pred, with, pX)
g = sapply(pred, with, pY)
r = f/g 
q_r = apply(r, 1, function(x) quantile(x, c(0.025, 0.975)))
y = c(q_r[1, ], rowMeans(r), q_r[2,])
z = rep(c("2.5%", "mean", "97.5%"), each = length(Grid))
df_r = data.frame(x = Grid,
                  y = y, 
                  z = z)
                  
df_r = df_r %>% filter(x > min(c(X,Y)), x < max(c(X,Y)))

# scale_y_log10 is ignoring my break at 0 for whatever reason
p1 = ggplot(df_r) +
  aes(x = x, y = y, linetype = z) +
  geom_line() +
  scale_y_log10() +
#  scale_y_log10(breaks = c(0, 1, 10, 100), 
#                labels  = c(0, 1, 10, 100)) +
  scale_linetype_manual(values = c("2.5%" = "dashed",
                                   "mean" = "solid",
                                   "97.5%" = "dashed")) +
  theme_tufte() +
  theme(legend.position = "none") +
  ylab("Density ratio (Control / Germ-free)") +
  xlab("Time of death (days)") 


#########
# f & g #
#########

# getting quantiles for intervals
q_f = apply(f, 1, function(x) quantile(x, c(0.025, 0.975)))
q_g = apply(g, 1, function(x) quantile(x, c(0.025, 0.975)))

# combining with posterior means
y = c(q_f[1, ], rowMeans(f), q_f[2,], 
      q_g[1, ], rowMeans(g), q_g[2,])

# formatting data.frame
z = rep(rep(c("2.5%", "mean", "97.5%"), each = length(Grid)), 2)
dens = rep(c("control", "germ-free"), each = 3*length(Grid))
df_fg = data.frame(x = c(Grid, Grid),
                   y = y, 
                   z = z,
                   dens = dens)

df_fg = df_fg %>% filter(x > min(c(X,Y)), x < max(c(X,Y)))

df_cont = df_fg %>% filter(dens == "control")
df_germ = df_fg %>% filter(dens != "control")

df_x = as.data.frame(X)
df_y = as.data.frame(Y) 

ylims = c(0, 0.003)
pf = ggplot() +
  geom_line(data = df_cont, aes(x = x, y = y, linetype = z)) +
  geom_histogram(data = df_x, aes(x = X,  y = ..density..), alpha = 0.2) +
  scale_linetype_manual(values = c("2.5%" = "dashed",
                                   "mean" = "solid",
                                   "97.5%" = "dashed")) +
  theme_tufte() +
  theme(legend.position = "none") +
  ylim(ylims) +
  xlab("Time of death (days)") +
  ylab("Density (control)")


pg = ggplot() +
  geom_line(data = df_germ, aes(x = x, y = y, linetype = z)) +
  geom_histogram(data = df_y, aes(x = Y,  y = ..density..), alpha = 0.2) +
  scale_linetype_manual(values = c("2.5%" = "dashed",
                                   "mean" = "solid",
                                   "97.5%" = "dashed")) +
  theme_tufte() +
  theme(legend.position = "none") +
  ylim(ylims) +
  xlab("Time of death (days)") +
  ylab("Density (germ-free)")


# arranging figure
grid.arrange(p1, pf, pg, layout_matrix = matrix(c(1, 2,
                                                  1, 3), nrow = 2, byrow = 2))


