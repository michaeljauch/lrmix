## Function to compute the stick-breaking
sb =  function(v) exp(c(log(v),0)+cumsum(log(c(1,1-v))))

#### The bnp_l1 runs the mcmc described in Section 3.4 of the paper and 
#### Section S2 of the supplementary material
# Input
# prior: a list of values for the hyperparameter of the model
# mcmcm: a list of values defining the number of mcmc iterations
# Grid:  grid of values for density estimation
## Output
# mcmc.chain: nsave draws (after burning and thinning) from the posterior 
#             distribution of (w1, v1, mu1, s21, w2, v2, mu2, s22, \theta)
# dens:       nsave draws (after burning and thinning) from the posterior 
#             distribution of (f,g,\theta) 
bnp_lr <- function(X, Y, prior, mcmc.par, Grid)
{
  # Hyper-parameters
  L = prior$L
  
  alpha = prior$alpha
  
  mU = prior$m
  lambdaU = prior$c
  aU = prior$a1
  bU = prior$a2
  
  mY = prior$m
  lambdaY = prior$c
  aY = prior$a1
  bY = prior$a2
  
  a1 = prior$b1
  a2 = prior$b2
  Kappa = prior$p0
  
  ## Initial values
  Q = 1
  
  R = sample(c(1,1),M,replace = T)
  
  muY = rnorm(L,0,1)
  s2Y = runif(L,0.5,2)
  vY = rbeta(L-1,1,alpha)
  wY = sb(vY)
  
  muX = rnorm(L,0,1)
  s2X = runif(L,0.5,2)
  vX = rbeta(L-1,1,alpha)
  wX = sb(vX)
  
  U = 1.001*X
  muU = rnorm(L,0,1)
  s2U = runif(L,0.5,2)
  vU = rbeta(L-1,1,alpha)
  wU = sb(vU)
  
  theta = 0.5
  
  Jchain = 0
  ############################################################
  ## By running an MCMC this script find "good" initial values
  ## for vY, muY, and s2Y.
  ############################################################
  ###########################################
  # Initial values for fY
  ###########################################
  # pb = txtProgressBar(min = 0, max = 10000, initial = 0, style = 4)
  cat("\n")
  cat("Finding initial g", sep="\n")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = 10000, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar  # Run MCMC
  out  =  list()
  z = 0
  for(s in 1:10000)
  {
    # Update vY muY s2Y
    logpost = function(x0,Y)
    {
      vY = x0[1:(L-1)]
      muY = x0[L:(2*L-1)]
      s2Y = x0[(2*L):(3*L-1)]
      
      wY = sb(vY)
      tmpfY = sum(log(sapply(1:L,function(l)dnorm(Y,muY[l],sqrt(s2Y[l])))%*%wY))
      tmpPrior = sum(sapply(1:L,function(l)dnorm(muY[l],mY,sqrt(s2Y/lambdaY)[l],log = T))) + 
        sum(LaplacesDemon::dinvgamma(s2Y,aY,bY,log = T))
      
      tmpfY + tmpPrior
    }
    
    x1 = x0 = c(vY,muY,s2Y)
    w = c(rep(0.1,L-1),rep(0.5,L),rep(0.5,L))
    
    ## Slice sampler
    # Step a)
    fx0 = log(runif(1))+logpost(x0,Y)
    
    # Step b)
    u = runif(3*L-1)
    LL = x0-w*u
    LL[1:(L-1)] = ifelse(LL[1:(L-1)]<0,0,LL[1:(L-1)])
    LL[(2*L):(3*L-1)] = ifelse(LL[(2*L):(3*L-1)]<0,0,LL[(2*L):(3*L-1)])
    
    RR = LL+w
    RR[1:(L-1)] = ifelse(RR[1:(L-1)]>1,1,RR[1:(L-1)])
    
    # Step c)
    x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    
    while(fx0 >= logpost(x1,Y))
    {
      LL = ifelse(x1<x0,x1,LL)
      RR = ifelse(x1>=x0,x1,RR)
      x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    }
    
    vY = x1[1:(L-1)]
    muY = x1[L:(2*L-1)]
    s2Y = x1[(2*L):(3*L-1)]
    wY = sb(vY)
    
    if(s>2000 & (s-2000)%%8==0)
    {
      z = z+1
      out[[z]] = list(vY = vY, muY = muY, s2Y = s2Y, wY = wY)
    }
    
    # Print progress
    setTxtProgressBar(pb,s)
  }
  
  
  PredY = 
    sapply(1:length(out),function(z)
      sapply(1:L,function(l)
        dnorm(Grid,out[[z]]$muY[l],sqrt(out[[z]]$s2Y[l])))%*%out[[z]]$wY)
  
  EpredY = rowMeans(PredY)
  Ind = which.min(sapply(1:length(out),function(z) mean(abs(PredY[,z] - EpredY)^2)))
  
  vY = out[[Ind]]$vY
  muY = out[[Ind]]$muY
  s2Y = out[[Ind]]$s2Y
  wY = out[[Ind]]$wY
  predY = sapply(1:L,function(l)dnorm(Grid,muY[l],sqrt(s2Y[l])))%*%wY
  

  ############################################################
  ## By running an MCMC this script find a "good" initial values
  ## for theta, vU, muU, and s2U. Based on the MCMC, this script also find a
  ## the pseudo prior for theta, vU, muU, and s2U.
  ############################################################
  ###########################################
  # Pseudo-prior for fU and theta and initial values  
  ###########################################
  cat("\n \n")
  cat("Finding pseudo priors and initial u", sep="\n")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = 10000, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar  # Run MCMC
  
  # Run MCMC
  MU = c()
  Theta = c()
  Theta[1] = theta
  
  U = seq(-5,5,len = 1001)
  cte = sapply(1:L,function(l)pnorm(U,muY[l],sqrt(s2Y[l])))%*%wY
  dyx = sapply(1:M,function(j)
    as.vector(sapply(1:L,function(l) dnorm(X[j],muY[l],sqrt(s2Y[l])))%*%wY))
  z=0
  
  for(s in 1:10000)
  {
    # Update vU muU s2U
    logpost = function(x0,X)
    {
      vU = x0[1:(L-1)]
      muU = x0[L:(2*L-1)]
      s2U = x0[(2*L):(3*L-1)]
      
      wU = sb(vU)
      
      
      densU  = sapply(1:L,function(l)dnorm(U,muU[l],sqrt(s2U[l])))%*%wU
      
      logdensx = function(j)
      {
        dx = (X[j]<=U)*dyx[j]
        log(theta*dyx[j]+(1-theta)*sum(densU*dx*(U[2]-U[1])/cte))
      }
      # tmp1fX = simplify2array(mclapply(1:M,logdensx,mc.cores = 1))
      tmp1fX = simplify2array(lapply(1:M,logdensx))
      
      tmpPrior = sum(sapply(1:L,function(l)dnorm(muU[l],mU,sqrt(s2U/lambdaU)[l],log = T))) + 
        sum(LaplacesDemon::dinvgamma(s2U,aU,bU,log = T))
      
      sum(tmp1fX) + tmpPrior
    }
    
    x1 = x0 = c(vU,muU,s2U)
    w = c(rep(0.1,L-1),rep(0.5,L),rep(0.5,L))
    
    ## Slice sampler
    # Step a)
    fx0 = log(runif(1))+logpost(x0,X)
    
    # Step b)
    u = runif(3*L-1)
    LL = x0-w*u
    LL[1:(L-1)] = ifelse(LL[1:(L-1)]<0,0,LL[1:(L-1)])
    LL[(2*L):(3*L-1)] = ifelse(LL[(2*L):(3*L-1)]<0,0,LL[(2*L):(3*L-1)])
    
    RR = LL+w
    RR[1:(L-1)] = ifelse(RR[1:(L-1)]>1,1,RR[1:(L-1)])
    
    # Step c)
    x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    
    # J = 0
    while(fx0 >= logpost(x1,X))
    {
      # J=J+1
      # print(paste("J in slide = ",J))
      LL = ifelse(x1<x0,x1,LL)
      RR = ifelse(x1>=x0,x1,RR)
      x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    }
    
    vU = x1[1:(L-1)]
    muU = x1[L:(2*L-1)]
    s2U = x1[(2*L):(3*L-1)]
    wU = sb(vU)
    
    # Update R and theta
    # aux = t(sapply(1:M,function(j){
    #   c(log(theta)+dnorm(X[j],muY[K[j]],sqrt(s2Y[K[j]]),log=T),
    #     log(1-theta)+dnorm(X[j],muY[K[j]],sqrt(s2Y[K[j]]),log=T)-
    #       pnorm(U[j],muY[K[j]],sqrt(s2Y[K[j]]),log = T))}))
    # R = apply(aux,1,function(tmp) sample(1:0,1,prob=exp(tmp-max(tmp))/sum(exp(tmp-max(tmp)))))
    
    densU  = sapply(1:L,function(l)dnorm(U,muU[l],sqrt(s2U[l])))%*%wU
    
    logdensx = function(j)
    {
      dx = (X[j]<=U)*dyx[j]
      log(c(theta*dyx[j],(1-theta)*sum(densU*dx*(U[2]-U[1])/cte)))
    }
    # probR = t(simplify2array(mclapply(1:M,logdensx,mc.cores = 1)))
    probR = t(simplify2array(lapply(1:M,logdensx)))
    R = apply(probR,1,function(tmp) sample(1:0,1,prob=exp(tmp-max(tmp))/sum(exp(tmp-max(tmp)))))
    
    
    theta = rbeta(1, a1 + sum(R), a2+M-sum(R))
    # theta = 0.25
    
    # if(s%%1000==0 | s==1) print(paste("chain",Jchain,"--", s,"iteration of 10000 (finding init. values and ps-prior for theta and fU)"))  
    
    if(s>200 & (s-200)%%8==0)
    {
      z = z+1
      Theta[z] = theta
      out[[z]] = list(vU = vU, muU = muU, s2U = s2U, wU = wU)
    }
    # Print progress
    setTxtProgressBar(pb,s)
  }
  
  
  # pseudo prior for base measure and precision
  Ind = sapply(1:length(out),function(z)sample(1:L,1,prob=out[[z]]$wU))
  mu_logs2_U = t(sapply(1:length(out),function(z) c(out[[z]]$muU[Ind[z]],log(out[[z]]$s2U[Ind[z]]))))
  
  
  library(ks)
  pseudo_prior_mu_logs2_U <- kde(x=mu_logs2_U)
  
  pseudo_alphaU = mean(sapply(1:length(out),function(z)
    optimize(function(alpha)sum(dbeta(out[[z]]$vU,1,alpha,log = T)),c(0,1e4),maximum = TRUE)$max))
  
  # pseudo prior fot theta
  auxmean = mean(Theta)
  #auxvar = 1.2*var(Theta) # 20% larger?
  auxvar = var(Theta)
  nu = ((auxmean*(1-auxmean))/auxvar)-1
  pseudo_a = auxmean*nu
  pseudo_b = (1-auxmean)*nu
  

  ###################################################
  ###################################################
  # MCMC
  ###################################################
  ####################################################

  out = list() # MCMC sequences are store here
  
  nburn = mcmc.par$nburn
  nsave = mcmc.par$nsave
  nskip = mcmc.par$nskip
  
  thetatmp = theta
  z = 0

  cat("\n \n")
  cat("Running MCMC", sep="\n")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = (nburn+nsave*nskip), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar  # Run MCMC
  
  for(s in 1:(nburn+nsave*nskip))
  {
    
    # Update vY muY s2Y
    densU  = sapply(1:L,function(l)dnorm(U,muU[l],sqrt(s2U[l])))%*%wU
    
    logpost = function(x0,Y,X)
    {
      vY = x0[1:(L-1)]
      muY = x0[L:(2*L-1)]
      s2Y = x0[(2*L):(3*L-1)]
      
      wY = sb(vY)
      tmpfY = sum(log(sapply(1:L,function(l)dnorm(Y,muY[l],sqrt(s2Y[l])))%*%wY))
      
      dyx = sapply(1:M,function(j)
        as.vector(sapply(1:L,function(l) dnorm(X[j],muY[l],sqrt(s2Y[l])))%*%wY))
      
      cte = sapply(1:L,function(l)pnorm(U,muY[l],sqrt(s2Y[l])))%*%wY
      
      logdensx = function(j)
      {
        dx = (X[j]<=U)*dyx[j]
        log(theta*dyx[j]+(1-theta)*sum(densU*dx*(U[2]-U[1])/cte))
      }
      # if(Q == 0)
      #   tmp1fX = sum(simplify2array(mclapply(1:M,logdensx,mc.cores = 1)))
      if(Q == 0)
        tmp1fX = sum(simplify2array(lapply(1:M,logdensx)))
      if(Q == 1)
        tmp1fX = sum(log(dyx))
      
      tmpPrior = sum(sapply(1:L,function(l)dnorm(muY[l],mY,sqrt(s2Y/lambdaY)[l],log = T))) + 
        sum(LaplacesDemon::dinvgamma(s2Y,aY,bY,log = T))
      
      tmpfY + tmp1fX + tmpPrior
    }
    
    x1 = x0 = c(vY,muY,s2Y)
    w = c(rep(0.1,L-1),rep(0.5,L),rep(0.5,L))
    
    ## Slice sampler
    # Step a)
    fx0 = log(runif(1))+logpost(x0,Y,X)
    
    # Step b)
    u = runif(3*L-1)
    LL = x0-w*u
    LL[1:(L-1)] = ifelse(LL[1:(L-1)]<0,0,LL[1:(L-1)])
    LL[(2*L):(3*L-1)] = ifelse(LL[(2*L):(3*L-1)]<0,0,LL[(2*L):(3*L-1)])
    
    RR = LL+w
    RR[1:(L-1)] = ifelse(RR[1:(L-1)]>1,1,RR[1:(L-1)])
    
    # Step c)
    x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    
    while(fx0 >= logpost(x1,Y,X))
    {
      LL = ifelse(x1<x0,x1,LL)
      RR = ifelse(x1>=x0,x1,RR)
      x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
    }
    
    vY = x1[1:(L-1)]
    muY = x1[L:(2*L-1)]
    s2Y = x1[(2*L):(3*L-1)]
    wY = sb(vY)
    
    # Update vU muU s2U
    if(Q == 0)
    {
      U = seq(-5,5,len = 1001)
      cte = sapply(1:L,function(l)pnorm(U,muY[l],sqrt(s2Y[l])))%*%wY
      dyx = sapply(1:M,function(j)
        as.vector(sapply(1:L,function(l) dnorm(X[j],muY[l],sqrt(s2Y[l])))%*%wY))
      
      logpost = function(x0,X)
      {
        vU = x0[1:(L-1)]
        muU = x0[L:(2*L-1)]
        s2U = x0[(2*L):(3*L-1)]
        
        wU = sb(vU)
        
        
        densU  = sapply(1:L,function(l)dnorm(U,muU[l],sqrt(s2U[l])))%*%wU
        
        logdensx = function(j)
        {
          dx = (X[j]<=U)*dyx[j]
          log(theta*dyx[j]+(1-theta)*sum(densU*dx*(U[2]-U[1])/cte))
        }
        # tmp1fX = simplify2array(mclapply(1:M,logdensx,mc.cores = 1))
        tmp1fX = simplify2array(lapply(1:M,logdensx))
        
        tmpPrior = sum(sapply(1:L,function(l)dnorm(muU[l],mU,sqrt(s2U/lambdaU)[l],log = T))) + 
          sum(LaplacesDemon::dinvgamma(s2U,aU,bU,log = T))
        
        sum(tmp1fX) + tmpPrior
      }
      
      x1 = x0 = c(vU,muU,s2U)
      w = c(rep(0.1,L-1),rep(0.5,L),rep(0.5,L))
      
      ## Slice sampler
      # Step a)
      fx0 = log(runif(1))+logpost(x0,X)
      
      # Step b)
      u = runif(3*L-1)
      LL = x0-w*u
      LL[1:(L-1)] = ifelse(LL[1:(L-1)]<0,0,LL[1:(L-1)])
      LL[(2*L):(3*L-1)] = ifelse(LL[(2*L):(3*L-1)]<0,0,LL[(2*L):(3*L-1)])
      
      RR = LL+w
      RR[1:(L-1)] = ifelse(RR[1:(L-1)]>1,1,RR[1:(L-1)])
      
      # Step c)
      x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
      
      while(fx0 >= logpost(x1,X))
      {
        LL = ifelse(x1<x0,x1,LL)
        RR = ifelse(x1>=x0,x1,RR)
        x1 = sapply(1:length(x0),function(i1)runif(1,LL[i1],RR[i1]))
      }
      
      vU = x1[1:(L-1)]
      muU = x1[L:(2*L-1)]
      s2U = x1[(2*L):(3*L-1)]
      wU = sb(vU)
    }
    if(Q == 1)
    {
      vU = rbeta(L-1,1,pseudo_alphaU)
      wU = sb(vU)
      aux = rkde(L, pseudo_prior_mu_logs2_U)
      muU = aux[,1]
      s2U = exp(aux[,2])
    }
    
    # Update R theta
    dyx = sapply(1:M,function(j)
      as.vector(sapply(1:L,function(l) dnorm(X[j],muY[l],sqrt(s2Y[l])))%*%wY))
    densU  = sapply(1:L,function(l)dnorm(U,muU[l],sqrt(s2U[l])))%*%wU
    cte = sapply(1:L,function(l)pnorm(U,muY[l],sqrt(s2Y[l])))%*%wY
    
    logdensx = function(j)
    {
      dx = (X[j]<=U)*dyx[j]
      log(c(theta*dyx[j],(1-theta)*sum(densU*dx*(U[2]-U[1])/cte)))
    }
    # probR = t(simplify2array(mclapply(1:M,logdensx,mc.cores = 1)))
    probR = t(simplify2array(lapply(1:M,logdensx)))
    
    if(Q==0) R = apply(probR,1,function(tmp) sample(1:0,1,prob=exp(tmp-max(tmp))/sum(exp(tmp-max(tmp)))))
    if(Q==1) R = rep(1,M)
    
    # Q
    theta = thetatmp
    logdensx = function(j)
    {
      dx = (X[j]<=U)*dyx[j]
      log(c(dyx[j],theta*dyx[j]+(1-theta)*sum(densU*dx*(U[2]-U[1])/cte)))
    }
    # aux0 = colSums(t(simplify2array(mclapply(1:M,logdensx,mc.cores = 1))))
    aux0 = colSums(t(simplify2array(lapply(1:M,logdensx))))
    aux = c(log(Kappa)+dbeta(theta,pseudo_a,pseudo_b, log = T), log(1-Kappa)+dbeta(theta, a1, a2, log = T))
    aux = ifelse(aux==Inf,1e6,aux)+aux0
    Q = rbinom(1,1,(exp(aux-max(aux))/sum(exp(aux-max(aux))))[1])
    
    # theta
    if(Q == 1) theta = rbeta(1, pseudo_a,pseudo_b)
    if(Q == 0) theta = rbeta(1, a1 + sum(R), a2+M-sum(R))
    
    thetatmp = theta
    theta = (1-Q)*thetatmp + Q
    # print(cbind(s,thetatmp,Q,theta))
    
    # if(s==1) print(paste("chain",Jchain,"--", "MCMC starts"))
    # if(s==1 | (s<=nburn & s%%100==0)) print(paste("chain",Jchain,"--", s,"iteration of", nburn+nsave*nskip, "Burning"))  
    # if(s>nburn & s%%100==0) print(paste("chain",Jchain,"--", s,"iteration of", nburn+nsave*nskip))  
    
    if(s %in% (nburn+(1:nsave)*nskip))
    {
      z = z + 1
      out[[z]] = list(theta=theta, 
                      vY = vY, muY = muY, s2Y = s2Y, wY = wY,
                      vU = vU, muU = muU, s2U = s2U, wU = wU)
    }
    # if(s>nburn & s%%100==0) save(out,file = paste("mcmc_seq/chain",Jchain,".RData",sep=""))
    # Print progress
    setTxtProgressBar(pb,s)
    
  }
  
  #######################
  # MCMCM post-processing: computing densities
  #######################
  # Compute fX and fX
  cat("\n \n")
  cat("MCMC post-processing: computing densities", sep="\n")
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(out), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar  # Run MCMC
  
  pred = list()
  out2 = list()
  
  # Function to compute fX and fY
  for(iter in 1:length(out))
  {
    theta = out[[iter]]$theta
    vY = out[[iter]]$vY
    muY = out[[iter]]$muY
    s2Y = out[[iter]]$s2Y
    wY = out[[iter]]$wY
    vU = out[[iter]]$vU
    muU = out[[iter]]$muU
    s2U = out[[iter]]$s2U
    wU = out[[iter]]$wU
    
    L = length(wU)
    
    predY = sapply(1:L,function(l)dnorm(Grid,muY[l],sqrt(s2Y[l])))%*%wY
    
    u = seq(-5,5,len = 10001)
    densU  = sapply(1:L,function(l)dnorm(u,muU[l],sqrt(s2U[l])))%*%wU
    cte = sapply(1:L,function(l)pnorm(u,muY[l],sqrt(s2Y[l])))%*%wY
    
    densx = function(j)
    {
      dx = (Grid[j]<=u)*as.vector(sapply(1:L,function(l) dnorm(Grid[j],muY[l],sqrt(s2Y[l])))%*%wY)
      dyx = as.vector(sapply(1:L,function(l) dnorm(Grid[j],muY[l],sqrt(s2Y[l])))%*%wY)
      theta*dyx+(1-theta)*sum(densU*dx*(u[2]-u[1])/cte)
    }
    tmp2fX <- simplify2array(lapply(1:length(Grid),densx))
    
    pY = predY#*(1/(Grid0*sd(Z)))
    pX = tmp2fX#*(1/(Grid0*sd(Z)))
    
    # plogY = predY*(1/(sd(Z)))
    # plogX = tmp2fX*(1/(sd(Z)))
    
    # if(iter==1 | iter%%100==0) print(paste(iter,"iteration of", length(out)))  
    # Print progress
    setTxtProgressBar(pb,iter)
    
    # list(pY = pY, plogY = plogY, pX = pX, plogX = plogX, plogY, theta = theta)
    out2[[iter]] = list(w1 = wY, v1 = vY, mu1 = muY, s21 = s2Y,
                        w2 = wU, v2 = vU, mu2 = muU, s22 = s2U,
                        theta = theta)
    pred[[iter]] = list(f = pX, g = pY, theta = theta)
  }
  cat("\n")

  list(mcmc.chain = out2, dens = pred)
}
