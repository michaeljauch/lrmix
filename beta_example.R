library(tidyverse)
library(ggthemes)
library(gridExtra)
# fx = Beta(a,b), fy = Beta(a, b-delta)
a = 1
b = 3
delta = 2
u = seq(0.01, 0.99, 0.001)
theta = beta(a,b-delta)/beta(a,b)
neg_dr = -delta*(1-u)^(-1+delta)*beta(a,b-delta)/beta(a,b)
fu = pbeta(u,a, b-delta)*neg_dr/(1-theta)
fx = dbeta(u, a, b)
fy = dbeta(u, a, b-delta)
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
