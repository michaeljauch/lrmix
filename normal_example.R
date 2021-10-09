library(tidyverse)
library(ggthemes)
library(gridExtra)
mu = -1
delta = 2
sigma = 1
u = seq(-5, 5, 0.01)
fu = pnorm((u-mu-delta)/sigma)*delta/(sigma^2)*exp(delta*(delta - 2*(u-mu))/(2*sigma^2))
fx = dnorm(u, mean = mu, sd = sigma)
fy = dnorm(u, mean = mu+delta, sd = sigma)
df = data.frame(f = c(fx, fy),
                x = rep(u, 2), 
                dist = rep(c("F","G"), each = length(u)))
  

p1 = ggplot(df) +
  aes(x = x, y = f, color = dist, linetype = dist) +
    geom_line(size = 1.2) +
    ylab("Density") +
    scale_color_manual(values = c("azure3", "azure4")) +
    scale_linetype_manual(values = c("solid", "solid")) +
    theme_tufte() +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(expression("f and g"))

df1 = data.frame(f = fu,
                 x = u,
                 dist = rep("U", length(u)))

p2 = ggplot(df1) +
  aes(x = x, y = f, color = dist, linetype = dist) +
  geom_line(size = 1.2) +
  ylab("Density") +
  scale_color_manual(values = "black") +
  scale_linetype_manual(values = "solid") +
  theme_tufte() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("u")
  

grid.arrange(p1, p2, nrow = 2)


