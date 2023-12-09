library(ggplot2)
library(tidyverse)
library(ggthemes)
library(gridExtra)
# F = Beta(1, 3) and G = Beta(1, 1)
u = seq(0.01, 0.99, 0.001)

fu = dbeta(u, 2 ,2)
fx = dbeta(u, 1, 3)
fy = dbeta(u, 1, 1)
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
