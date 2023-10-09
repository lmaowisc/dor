

# ICM for monotonization
CM=function(G,Q){
  m=length(G)
  y=rep(NA,m-1)
  # Q=c(0,Q)
  # G=c(0,G)
  j=1
  while(j<m){
    slopes=rep(NA,m-j)
    for (k in (j+1):m){
      slopes[k-j]=(Q[k]-Q[j])/(G[k]-G[j])
    }
    j1=which.min(slopes)+j
    y[j:(j1-1)]=min(slopes)
    j=j1
  }
  return(y)
}

#######################################################
#   Main function: Naive and improved IPCW estimator  #
#######################################################
# Requires: stats 
# Input: (x1, delta1, x2, delta2, tau)
#   x1: time to earliest of response, outcome event, or censoring 
#   delta1: indicator of response (or of response/event)
#   x2: time to earlier of outcome event and censoring 
#   delta2: indicator of outcome event
#   tau: restriction time
#   mono: whether to monotonized improved IPW
#   M: # number of time knots to evaluate improved IPW 
#
# Output: (t1,surv1, se1, t2, surv2, surv2m, se2)
#   t1,surv1, se1: times, naive estimates, standard errors
#   t2,surv2, se2: times, improved IPW estimates, standard errors
#   surv2m: monotonized improved IPW (if mono=T)
#############################################################
dorfit <- function(x1, delta1, x2, delta2, tau, mono=TRUE, M=100)
{
  
  ##############################################################
  #            Preparations                                    #
  ##############################################################
  # sample size 
  n <- length(x1)
  
  # refined the data to incorporate the
  # restriction time tau
  x1[x1>x2] <- x2[x1>x2]
  x1 <- pmin(x1, tau)
  delta1[x1==tau] <- 1
  x2 <- pmin(x2, tau)
  delta2[x2==tau] <- 1
  
  ## define delta1 to indicate response or event without response
  delta1[delta2==1] <- 1
  
  # estimate censoring survival function
  # G_C(t)
  deltaC <- 1-delta2 # indicator for censoring
  # fit KM curves for G_C(t)
  fitC <- survfit(Surv(x2, deltaC) ~ 1)
  tC <- c(0, fitC$time)
  survC <- c(1, fitC$surv)
  cdfC <- ecdf(tC) # empr dist of timeC 
  # can be used to get G_C(t) from survC for any t
  nC <- length(tC)
  if (survC[nC]==0){survC[nC]=survC[nC-1]}
  
  
  Dt <- pmax(x2 - x1, 0) # should require x>=x1 to begin with
  # ordered unique values of Dt
  t1 <- sort(unique(c(0,Dt)))
  m1 <- length(t1)
  
  ##############################################################
  #            Naive IPCW                                      #
  ##############################################################
  
  ## evaluate G_C(x2) ##
  x2r <- round(cdfC(x2)*nC) # ranks of x2 in tC
  GCx2 <- survC[x2r] # n-vector of G_C(x2) 
  
  DtMAT <- outer(Dt,t1,">") # n by m matrix of I(Dt > t1)
  # dim(DtMAT)
  # ep <- 1e-12
  
  # naive IPW
  surv1MAT <- matrix(delta2/GCx2, n, m1)*DtMAT
  surv1 <- colMeans(surv1MAT)
  
  
  ##############################################################
  #            Improved IPCW                                   #
  ##############################################################
  
  # Vector of t2 at which S_IPW changes value
  tdn <- unique(pmax(as.vector(outer(tC,x1[delta1==1],"-")),0))
  ## t's at which G_C(X_1+t) changes value
  
  ## combine with t1
  t2 <- unique(sort(c(t1,tdn)))
  m2 <- length(t2)
  
  if (m2>M){
    t2 <- seq(0,tau,length=M)
    m2 <- M
  }
  
  
  
  x1tMAT <- outer(x1,t2,"+")
  
  x1tMATr <- round(cdfC(x1tMAT)*nC)
  
  GCx1t <- matrix(survC[x1tMATr],n,m2)
  
  
  DtMAT2 <- outer(Dt,t2,">") # n by m matrix of I(Dt > t1)
  surv2MAT <- DtMAT2/GCx1t
  
  surv2 <- colMeans(surv2MAT)
  
  if (mono){
    t2m <- c(-0.1,t2)
    
    Q <- c(0,cumsum((1-surv2)*diff(t2m)))
    surv2m <-  1- CM(t2m,Q)
    
    surv2m <- pmin(1,pmax(surv2m,0))
  }else{
    surv2m <- NULL
  }
  ##############################################################
  ###    Variance estimation for naive/improved IPCW         ###
  ##############################################################
  
  
  ## censoring martingale M_C
  Ut <- sort(unique(x2[deltaC==1]))
  l <- length(Ut)
  
  eeq <- function(y1,y2){
    abs(y1-y2)<1e-8
  }
  
  dNC <- outer(x2,Ut,eeq)*matrix(rep(deltaC,l),n,l)
  
  R2MAT <- outer(x2,Ut,">=")
  pi2 <- colMeans(R2MAT)
  
  dLambdaC <- colMeans(dNC)/pi2
  
  dMC <- dNC - R2MAT*matrix(rep(dLambdaC,each=n),n,l)
  
  # colMeans(dMC)
  ## matrices hold the IF's of surv1 and surv2
  IF1CMAT <- matrix(0,n,m1)
  IF2CMAT <- matrix(0,n,m2)
  
  
  # j <- 1
  for (j in 1:l){
    
    ## naive estimator
    zetatu <-  colMeans(surv1MAT*matrix(rep(R2MAT[,j],m1),n,m1))
    int1tu <- zetatu/pi2[j]
    IF1CMAT <- IF1CMAT + outer(dMC[,j],int1tu)
    
    ## improved IPW
    kappatu <- colMeans(surv2MAT*(x1tMAT>=Ut[j]))
    int2tu <- kappatu/pi2[j]
    IF2CMAT <- IF2CMAT + outer(dMC[,j],int2tu)
  }
  
  IF1MAT <- surv1MAT - matrix(rep(surv1,each=n),n,m1) + IF1CMAT
  IF2MAT <- surv2MAT - matrix(rep(surv2,each=n),n,m2) + IF2CMAT
  
  
  V1 <- colMeans(IF1MAT^2)/n
  se1 <- sqrt(V1)
  
  V2 <- colMeans(IF2MAT^2)/n
  se2 <- sqrt(V2)
  
  # sqrt(V1)
  
  
  
  
  
  
  
  return(list(t1=t1, surv1=surv1,t2=t2,surv2=surv2,surv2m=surv2m, se1=se1,se2=se2))
}







# Simulate uncensored bivariate response and event times
# from a log-normal distribution with normal mean mu and variance Sigma
# multiplied by a bivariate factor of gamma
bivsimu <- function(n, mu, Sigma,gamma){
  
  Tt <- mvrnorm(n, mu, Sigma)
  
  T1 <- exp(Tt[,1])*gamma[1]
  T2 <- exp(Tt[,2])*gamma[2]
  
  T1[T1>T2] <- T2[T1>T2]
  
  return(list(T1=T1,T2=T2))
  
}


### monte carlo to get the true S_D(t)
trueSD <- function(Ns, n0, mu, Sigma, gamma, tau, filename, M=500){
  
  t0 <- seq(0,tau,length=M)
  
  tmp <- matrix(NA,Ns,M)
  
  for (j in 1:Ns){
    dat0 <- bivsimu(n=N0, mu, Sigma,gamma)
    D0 <- pmin(dat0$T2,tau) - pmin(dat0$T1,tau)
    
    tmp[j,] <- 1 - ecdf(D0)(t0)
  }
  
  results <- matrix(NA,2,M)
  results[1,] <- t0
  results[2,] <- colMeans(tmp)
  
  write.table(results,filename)
}