library(progress)

set.seed(123)

## Read dataset
d = read.csv("Hoel1972.txt")
dim(d)

## Data in orignal scale (0,\infty)
Y0 = d[d[,1]=="Germ-free",3]
X0 = d[d[,1]=="Control",3]

## Transformation using the log
Y  = log(Y0)
X  = log(X0)

N = length(Y)
M = length(X)

## Standardize the transform data (This helps to specify default prior distributions) 
Z = c(Y,X)
Y = (Y - mean(Z))/sd(Z)
X = (X - mean(Z))/sd(Z)

# grid in the original and transform scale
Grid = seq(-5,6,len=1501)
Grid0 = exp(sd(Z)*Grid + mean(Z))
Grid = Grid[Grid0<1500]
Grid0 = Grid0[Grid0<1500]

## Prior
prior = list(L = 25, # Truncation level for stick-breaking
             alpha = 1, # precision parameter
             m = 0, c = 1, a1 = 1, a2 = 1, # base measure: normal-inv gamma(m,c,a1,a2)
             b1 = 1, b2 = 1, # \tilde{\theta} ~ beta(b1, b2) 
             p0 = 0.5) # \gamma ~ bernoulli(p0)

## mcmc -- number of iterations =  nburn + nskip*nsave 
mcmc.par = list(nburn = 2000, # number of iterations to be burned
                nskip = 10, # thinning parameter
                nsave = 1000) # number of iterations to be saved


## Fitting the model and running mcmc
source("MCMCalgorithm.R")
out = bnp_lr(X, Y, prior, mcmc.par, Grid)

## Summarizing the mcmc results
PDF = TRUE
if(PDF==TRUE)  pdf("Figures.pdf")

dens.f = sapply(out$dens, function(d) d$f)*(1/(sd(Z)))
dens.g = sapply(out$dens, function(d) d$g)*(1/(sd(Z)))
theta = sapply(out$dens, function(d) d$theta)


# Results for f -- log scale
dens.f.m = apply(dens.f,1,mean)
dens.f.l = apply(dens.f,1,quantile,prob=0.05)
dens.f.u = apply(dens.f,1,quantile,prob=0.95)
plot(x = log(Grid0), y = dens.f.m, type = "l", ylim = range(dens.f,dens.g), xlab = "log(x)", ylab = "f")
lines(x = log(Grid0), y = dens.f.l, lty = 2)
lines(x = log(Grid0), y = dens.f.u, lty = 2)
hist(log(X0), add = T, border = "gray", freq = F, breaks = 12)


# Results for g -- log scale
dens.g.m = apply(dens.g,1,mean)
dens.g.l = apply(dens.g,1,quantile,prob=0.05)
dens.g.u = apply(dens.g,1,quantile,prob=0.95)
plot(x = log(Grid0), y = dens.g.m, type = "l", ylim = range(dens.f,dens.g), xlab = "log(x)", ylab = "g")
lines(x = log(Grid0), y = dens.g.l, lty = 2)
lines(x = log(Grid0), y = dens.g.u, lty = 2)
hist(log(Y0), add = T, border = "gray", freq = F, breaks = 12)

# Results for r(x) =  g(x)/f(x) -- log scale
r = dens.g/dens.f
r.m = apply(r,1,mean)
r.l = apply(r,1,quantile,prob=0.05)
r.u = apply(r,1,quantile,prob=0.95)
plot(x = log(Grid0), y = r.m, type = "l", ylim = c(0,40), xlab = "log(x)", ylab = "r")
lines(x = log(Grid0), y = r.l, lty = 2)
lines(x = log(Grid0), y = r.u, lty = 2)

# Results for theta
hist(theta, xlim = c(0,1), freq=F, xlab = "theta", main = paste("P[theta = 1 | Data] = P[f = g | Data] =",round(mean(theta==1))))

if(PDF==TRUE)  dev.off()
