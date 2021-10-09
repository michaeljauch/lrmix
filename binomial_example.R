logit = function(p) {
  log(p/(1-p))
}

inv_logit = function(theta) {
  exp(theta)/(1+exp(theta))
}

pmf_u = function(u, n, theta, phi, delta) {
    pbinom(u, n, inv_logit(phi+delta))*(exp(delta)-1)/(1-theta)*((1+exp(phi+delta))/(1+exp(phi)))^n*exp(-delta*(1+u))
}


n = 10
phi = logit(0.4)
delta = 1
theta = exp(-n*delta)*((1+exp(phi))/(1+exp(phi+delta)))^(-n)
fx = dbinom(0:n, n, inv_logit(phi))
fy = dbinom(0:n, n, inv_logit(phi+delta))
fu = c(pmf_u(0:(n-1), n, theta, phi, delta), 0)
sum(fu) # checking this adds up to 1

df = data.frame(f = c(fx, fy, fu),
                x  = rep(0:n, 3),
                pmf = factor(rep(c("F","G","U"), each = (n+1))))

df$pmf = factor(df$pmf, levels = c("F", "G", "U"))

library(tidyverse)
library(ggthemes)
ggplot(df) +
  aes(x = x, y = f, fill = pmf) +
  geom_col(position = "dodge") +
  facet_wrap(~ pmf, nrow = 3) +
  ylab("probability") +
  scale_fill_manual(values = c("azure3", "azure4", "aquamarine3")) +
  scale_x_continuous(breaks = 0:n) +
  theme_tufte() +
  theme(legend.position = "none") 



df2 = data.frame(f = c(fx, fy),
                x  = rep(0:n, 2),
                pmf = factor(rep(c("F","G"), each = (n+1))))

df2$pmf = factor(df2$pmf, levels = c("F", "G"))

p1 = ggplot(df2) +
  aes(x = x, y = f, fill = pmf) +
  geom_col(position = "dodge", width = 0.3) +
  ylab("Probability") +
  xlab("") +
  scale_x_continuous(breaks = 0:n) +
  ylim(c(0,0.28)) +
  theme_tufte() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("azure3", "azure4")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression("f and g"))

df3 =  data.frame(f = fu,
                  x  = 0:n, 
                  pmf = rep("U", length(fu)))

p2 = ggplot(df3) +
  aes(x = x, y = f, fill = pmf) +
  geom_col(position = "dodge", width = 0.15) +
  ylab("Probability") +
  xlab("") +
  scale_x_continuous(breaks = 0:n) +
  ylim(c(0,0.28)) +
  theme_tufte() +
  theme(legend.position = "none") +
  scale_fill_manual(values = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("u")

library(gridExtra)
grid.arrange(p1, p2, nrow =2)


# ggplot(df) +
#  aes(x = x, y = f, fill = pmf) +
#  geom_col(position = "dodge", width = 0.5) +
#  ylab("probability") +
#  scale_fill_manual(values = c("azure3", "azure4", "aquamarine3")) +
#  scale_x_continuous(breaks = 0:n) +
#  theme_tufte() # +
#  theme(legend.position = "none") 




