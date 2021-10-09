library(Ternary)
# https://cran.r-project.org/web/packages/Ternary/vignettes/Ternary.html

pdf(file='~/ternary_plot.pdf',
    width=10, height=10)

par(mar = rep(0, 4))

TernaryPlot(alab = expression('Mass at ' ~ x[1]), 
            blab = expression('Mass at ' ~ x[2]), 
            clab = expression('Mass at ' ~ x[3]), 
            grid.lines = 5, grid.lty = 'dotted',
            grid.minor.lines = 1, grid.minor.lty = 'dotted',
            axis.labels=seq(0,1,.2), lab.offset = .15,
            lab.cex=2, padding=.12, ticks.length = .015, 
            axis.cex=1.3)

middle_triangle <- matrix(c(
  .2, .2, .6,
  .2/.4, .2/.4, 0,
  1, 0, 0
), ncol = 3, byrow = TRUE)
TernaryPolygon(middle_triangle, col = 'grey', border = 'black')

# Add data points
data_points <- list(c(.2, .2, .6), 
                    c(.2/.4, .2/.4, 0), 
                    c(1, 0, 0)
)
AddToTernary(points, data_points, pch = 21, cex = 7, bg='white')

AddToTernary(text, data_points, 
             c(expression(g), expression(g^x[2]), expression(g^x[1])), cex = 1.6, font = 2)

dev.off()



