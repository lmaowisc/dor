library(survival)
library(MASS)
rm(list=ls())


setwd("C:/Users/lmao/Dropbox/Lu Mao Team Folder/Collab projects/LJ Wei/DOR")
source("DORfunctions.R")


# gamma <- c(1,1)

Sigma <- matrix(c(1.01,-0.12,-0.12,1.41),2,2)
mu <- c(0.91, 2.38)
tau <- 30


# filename <- "./simuDat/a.txt"


################################################################################
#     Skip the following code if true values are generated                     #
################################################################################
Ns <- 1000
n0 <- 10000
# (a) gamma = (1,1)
trueSD(Ns, n0, mu, Sigma, gamma=c(1,1), tau=30, filename = "./simuDat/a.txt")
# (b) gamma = (1,4)
trueSD(Ns, n0, mu, Sigma, gamma=c(1,2), tau=30, filename = "./simuDat/b.txt")
# (c) gamma = (4,1)
trueSD(Ns, n0, mu, Sigma, gamma=c(2,2), tau=30, filename = "./simuDat/c.txt")
# (d) gamma = (4,4)
trueSD(Ns, n0, mu, Sigma, gamma=c(4,4), tau=30, filename = "./simuDat/d.txt")
################################################################################






## LT's example ###
#(1,1) (1, 4)  (4, 1) (4, 4).
#(a)#
Gamma <- matrix(NA,4,2)
Gamma[1,] <- c(1,1)
Gamma[2,] <- c(1,2)
Gamma[3,] <- c(2,2)
Gamma[4,] <- c(4,4)

rownames(Gamma) <- c("a","b","c","d")
N <- 5000
tau <- 30
M <- 100
t0 <- seq(0,tau,length=M)
##############################################################################
#   Simulation begins - save results to paste0("./simuDat/",sn,"n",n,".txt") #
##############################################################################


tk <- c(0,24)

for (n in c(50,200,500,1000)){



# g <- 1


for (g in 1:4){
  
  sn <- rownames(Gamma)[g]
  
truefile <- paste0("./simuDat/",sn,".txt")
truedat <- read.table(truefile)

S0 <- unlist(truedat[2,round(ecdf(truedat[1,])(t0)*ncol(truedat))])

gamma <- Gamma[g,]

{
  
  results <- matrix(NA,12,M)
  
  SURV1 <- matrix(NA,N,M)
  SE1 <- matrix(NA,N,M)
  SURV2 <- matrix(NA,N,M)
  SE2 <- matrix(NA,N,M)
  SURV2m <- matrix(NA,N,M)

  
  for (j in 1:N){
    

    dat <- bivsimu(n, mu, Sigma,gamma)
    T1 <- dat$T1
    T2 <- dat$T2
    
    
    # censoring distribution
    C <- pmin(runif(n, 0, 36),rexp(n,rate=0.003))
     # C <- runif(n, 0, 132)
    
    # C <- pmin(runif(n, 0, 45),rexp(n,rate=0.003))
    
    
    x1 <- pmin(T1, C)
    delta1 <- 1*(T1<=C) 
    
    x2 <- pmin(T2, C)
    delta2 <- 1*(T2<=C) 
    
    if (j==1){
    cat(mean(delta1),
    mean(delta2),
    mean(delta1*(1-delta2)),"\n")
    }
    
    fit <- dorfit(x1, delta1, x2, delta2, tau)
    
    t1 <- fit$t1
    surv1 <- fit$surv1
    se1 <- fit$se1
    
    t2 <- fit$t2
    surv2 <- fit$surv2
    surv2m <- fit$surv2m
    se2 <- fit$se2
    
    
    m1 <- round(ecdf(t1)(t0)*length(t1))
    m2 <- round(ecdf(t2)(t0)*length(t2))
    
    za <- qnorm(0.975)
    
    SURV1[j,] <- surv1[m1]
    SE1[j,] <- se1[m1]
    SURV2[j,] <- surv2[m2]
    SE2[j,] <- se2[m2]
    SURV2m[j,] <- surv2m[m2]
    

    ##
    if (j%%1000 == 0){
      for (k in 1:2){
      k0 <- ecdf(t0)(tk[k])*M
      Sk0 <- unlist(S0[k0])
      
      cat(n,j,tk[k],
          round(mean(SURV1[,k0] - Sk0, na.rm=T),4),
          round(sd(SURV1[,k0],na.rm=T),4),
          round(mean(SE1[,k0],na.rm=T),4),
          round(mean(abs(SURV1[,k0] - Sk0)<=za*SE1[,k0],na.rm=T),4),
          
          round(mean(SURV2[,k0] - Sk0,na.rm=T),4),
          round(sd(SURV2[,k0],na.rm=T),4),
          round(mean(SE2[,k0],na.rm=T),4),
          round(mean(abs(SURV2[,k0] - Sk0)<=za*SE2[,k0],na.rm=T),4),
          
          round(mean(SURV2m[,k0] - Sk0,na.rm=T),4),
          round(sd(SURV2m[,k0],na.rm=T),4),
          round(mean(SE2[,k0],na.rm=T),4),
          round(mean(abs(SURV2m[,k0] - Sk0)<=za*SE2[,k0],na.rm=T),4),
          "\n")
      }
    }
    
    
    
    
  }
  
  # naive estimator
  results[1,] <- colMeans(SURV1) - S0
  results[2,] <- apply(SURV1,2,sd)
  results[3,] <- colMeans(SE1)
  
  #ll-transformed CI
  llS0 <- log(-log(S0))
  llSE1 <- -SE1/(SURV1*log(SURV1))
  llSE1[SURV1==0] <- 0
  results[4,] <- colMeans(abs(log(-log(SURV1)) - matrix(rep(llS0,each=N),N,M))<=za*llSE1,na.rm=T)
  
  
  # IPW estimator
  results[5,] <- colMeans(SURV2) - S0
  results[6,] <- apply(SURV2,2,sd)
  results[7,] <- colMeans(SE2)
  llSE2 <- - SE2/(SURV2*log(SURV2))
  llSE2[SURV2==0] <- 0
  results[8,] <- colMeans(abs(log(-log(SURV2)) - matrix(rep(llS0,each=N),N,M))<=za*llSE2,na.rm=T)
  
  # mIPW estimator
  results[9,] <- colMeans(SURV2m) - S0
  results[10,] <- apply(SURV2m,2,sd)
  results[11,] <- results[7,]
  llSE2m <- -SE2/(SURV2m*log(SURV2m))
  llSE2m[SURV2m==0] <- 0
  results[12,] <- colMeans(abs(log(-log(SURV2m)) - matrix(rep(llS0,each=N),N,M))<=za*llSE2m,na.rm=T)
  
}

write.table(results,paste0("./simuDat/",sn,"n",n,".txt"))

}
}




#########################################################################
#             Graphics                                                  #
#########################################################################

######### average lines ####################################


par(mfrow=c(2,4))
# par(mar = c(5, 4.5, 0.5, 2))
par(oma=c(2.2,2.5,4,0))
par(mar=c(4.1, 4.5, 1.1, 1.1))
# c(5.1, 4.1, 4.1, 2.1). 


# red50 <- rgb(255, 0, 0, max = 255, alpha = 125, names = "red50")
# green50 <- rgb(t(col2rgb("green")), max = 255, alpha = 125, names = "green50") 
# purple50 <- rgb(t(col2rgb("blue")), max = 255, alpha = 125, names = "purple50")

for (n in c(50,200)){
for (g in 1:4){
  
sn <- rownames(Gamma)[g]

truedat <- read.table(paste0("./simuDat/",sn,".txt"))
results <- read.table(paste0("./simuDat/",sn,"n",n,".txt"))

# plot the true and estimated curves
S0 <- unlist(truedat[2,round(ecdf(truedat[1,])(t0)*ncol(truedat))])

surv1 <- results[1,] + S0
surv2 <- results[5,] + S0
surv2m <- results[9,] + S0

plot(t0,S0,type='l',xlim=c(0,tau),ylim=c(0,1),lwd=4,lty=3,
     # main=paste0("(",sn,")"),
     xlab="t",ylab=expression(S[D](t)),
     cex.lab=1.2,cex.axis=1.2,cex.main=1.2)
lines(t0,surv1,lwd=2,lty=1,col="red")
lines(t0,surv2,lwd=2,lty=1,col="blue")
lines(t0,surv2m,lwd=2,lty=1,col="purple")

# legend("topright",lty=c(3,1,1,1),col=c("black","red","blue","purple"),lwd=2,
#        c("True","Naive","IPW","IPW-mnt"),cex=1.09)

}
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',lty=c(3,1,1,1),col=c("black","red","blue","purple"),lwd=2,
       c("True","Initial","IPW","IPW-mnt"),cex=1.4,
       xpd = TRUE, horiz = TRUE,  seg.len=2, bty = 'n')

par(oma=c(2.2,2.5,4,0))
par(mar=c(4.1, 4.5, 1.1, 1.1))

mtext( '(a)', side=3, line=1, at=0.15, outer=TRUE,cex=1)
mtext( '(b)', side=3, line=1, at=0.4, outer=TRUE,cex=1)
mtext( '(c)', side=3, line=1, at=0.65, outer=TRUE,cex=1)
mtext( '(d)', side=3, line=1, at=0.9, outer=TRUE,cex=1)

mtext( 'n = 200', side=2, line=1, at=0.29, outer=TRUE,cex=1)
mtext( 'n = 50', side=2, line=1, at=0.79, outer=TRUE,cex=1)

###############################################################


###### plot some typical realizations #############

# n <- 50

par(mfrow=c(3,4))
# par(mar = c(5, 4.5, 0.5, 2))
par(oma=c(2.2,2.5,4,0))
par(mar=c(4.1, 4.5, 1.1, 1.1))
# c(5.1, 4.1, 4.1, 2.1). 
for (n in c(50,200,500)){
  for (g in 1:4){
    
    sn <- rownames(Gamma)[g]
  
    gamma <- Gamma[g,]
    truedat <- read.table(paste0("./simuDat/",sn,".txt"))
    # plot the true and estimated curves
    S0 <- unlist(truedat[2,round(ecdf(truedat[1,])(t0)*ncol(truedat))])
    
    dat <- bivsimu(n, mu, Sigma,gamma)
    T1 <- dat$T1
    T2 <- dat$T2
    
    
    # censoring distribution
    C <- pmin(runif(n, 0, 36),rexp(n,rate=0.003))
    # C <- runif(n, 0, 132)
    
    # C <- pmin(runif(n, 0, 45),rexp(n,rate=0.003))
    
    
    x1 <- pmin(T1, C)
    delta1 <- 1*(T1<=C) 
    
    x2 <- pmin(T2, C)
    delta2 <- 1*(T2<=C) 
    

    
    fit <- dorfit(x1, delta1, x2, delta2, tau)
    
    t1 <- fit$t1
    surv1 <- fit$surv1
    se1 <- fit$se1
    
    t2 <- fit$t2
    surv2 <- fit$surv2
    surv2m <- fit$surv2m
    se2 <- fit$se2
    
    
    m1 <- round(ecdf(t1)(t0)*length(t1))
    m2 <- round(ecdf(t2)(t0)*length(t2))
    

    plot(t0,S0,type='l',xlim=c(0,tau),ylim=c(0,1),lwd=4,lty=3,
         # main=paste0("(",sn,")"),
         xlab="t",ylab=expression(S[D](t)),
         cex.lab=1.2,cex.axis=1.2,cex.main=1.2)
    
    lines(t0,surv1[m1],lwd=2,lty=1,col="red")
    lines(t0,surv2[m2],lwd=2,lty=1,col="blue")
    lines(t0,surv2m[m2],lwd=2,lty=1,col="purple")
    
    # legend("topright",lty=c(3,1,1,1),col=c("black","red","blue","purple"),lwd=2,
    #        c("True","Naive","IPW","IPW-mnt"),cex=1.09)
    

  }
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',lty=c(3,1,1,1),col=c("black","red","blue","purple"),lwd=2,
       c("True","Initial","IPW","IPW-mnt"),cex=1.4,
       xpd = TRUE, horiz = TRUE,  seg.len=2, bty = 'n')
# xpd = TRUE makes the legend plot to the figure

par(oma=c(2.2,2.5,4,0))
par(mar=c(4.1, 4.5, 1.1, 1.1))

mtext( '(a)', side=3, line=1, at=0.15, outer=TRUE,cex=1)
mtext( '(b)', side=3, line=1, at=0.4, outer=TRUE,cex=1)
mtext( '(c)', side=3, line=1, at=0.65, outer=TRUE,cex=1)
mtext( '(d)', side=3, line=1, at=0.9, outer=TRUE,cex=1)

mtext( 'n = 500', side=2, line=1, at=0.195, outer=TRUE,cex=1)
mtext( 'n = 200', side=2, line=1, at=0.53, outer=TRUE,cex=1)
mtext( 'n = 50', side=2, line=1, at=0.86, outer=TRUE,cex=1)

################


###########################################################################
#             Tables for inferential results                              #
###########################################################################


##########################
# response and event rates
##########################

rates <- matrix(NA,4,2)
N <- 20000
for (g in 1:4){

  gamma <- Gamma[g,]


  dat <- bivsimu(N, mu, Sigma,gamma)
  T1 <- dat$T1
  T2 <- dat$T2
  C <- pmin(runif(N, 0, 36),rexp(n,rate=0.003))
  
  x1 <- pmin(T1, C)
  delta1 <- 1*(T1<=C) 
  
  x2 <- pmin(T2, C)
  delta2 <- 1*(T2<=C) 
  
  rates[g,1] <- mean(x1<x2)
  rates[g,2] <- mean(delta2)
}

n <- 500

tb <- c(0, 6, 12, 24)
mb <- length(tb)
nb <- 14

g <- 1


rdr <- function(x,r=3){format(round(x,r),nsmall=r)}
for (g in 1:4){
  
  
  sn <- rownames(Gamma)[g]
  gamma <- Gamma[g,]
  truedat <- read.table(paste0("./simuDat/",sn,".txt"))
  
  
  
  Sb <- unlist(truedat[2,round(ecdf(truedat[1,])(tb)*ncol(truedat))])
  
  # t0 <- unlist(truedat[1,])
  
  tmp <- read.table(paste0("./simuDat/",sn,"n",n,".txt"))
  
  dim(tmp)
  
  tmp1 <- t(tmp[,round(ecdf(t0)(tb)*length(t0))])
  
  tmp2 <- cbind(tmp1,(tmp1[,c(2,2)]/tmp1[,c(6,10)])^2)
  
  results <- cbind("&",rdr(tb,0),"&",rdr(Sb),"&&",
               rdr(tmp2[,1]), "&", rdr(tmp2[,2]), "&", rdr(tmp2[,3]), "&", rdr(tmp2[,4]), "&&",
               rdr(tmp2[,9]), "&", rdr(tmp2[,10]), "&", rdr(tmp2[,7]), "&", rdr(tmp2[,12]), "&&",
               rdr(tmp2[,14],2), "\\" 
  )
  
  
  
  
  rownames(results) <- rep("",mb)
  nb <- ncol(results)
  colnames(results) <- rep("",nb)
  results[mb,nb] <- "\\\noalign{\\medskip} "
  
  
  cat("\\multicolumn{14}{l}{",
    paste0("(",sn,")"), "Response rate:", paste0(round(100*rates[g,1]),"\\%; "),
      "Event rate:", paste0(round(100*rates[g,2]),"\\% } "),
    "\\\\")
  print(noquote(results))

  }

# (0.01909628/0.02041490)^2

# paste0("./simuDat/",sn,"n",n,".txt")












## simple example ###
{
  
  t0 <- 0
  n <- 500
  
  N <- 2000
  
  results <- matrix(NA,8,N)
  
  for (j in 1:N){
    c=runif(n, 0, 1.5)
    t1=rexp(n)
    t2=rexp(n)*1.5
    t1[t1>t2]=t2[t1>t2]
    M=100
    
    tau=1
    x1=pmin(t1, c)
    delta1=1*(t1<c)
    x2=pmin(t2, c)
    delta2=1*(t2<c)
    
    mean(delta1)
    mean(delta2)
    mean(delta1*(1-delta2))
    
    
    fit <- dorfit(x1, delta1, x2, delta2, tau)
    
    t1 <- fit$t1
    surv1 <- fit$surv1
    se1 <- fit$se1
    
    t2 <- fit$t2
    surv2 <- fit$surv2
    surv2m <- fit$surv2m
    se2 <- fit$se2
    
    
    m1 <- sum(t1<=t0)
    m2 <- sum(t2<=t0)
    
    za <- qnorm(0.975)
    
    results[1,j] <- surv1[m1]
    results[2,j] <- se1[m1]
    # results[3,j] <- abs()
    
    results[3,j] <- surv2[m2]
    results[4,j] <- se2[m2]
    
    cat(j,mean(results[1,],na.rm=T),sd(results[1,],na.rm=T),mean(results[2,],na.rm=T),
        mean(results[3,],na.rm=T),sd(results[3,],na.rm=T),mean(results[4,],na.rm=T),"\n")
    
  }
  
  
}




## example
############
n=50
c=runif(n, 0, 1.5)
t1=rexp(n)
t2=rexp(n)*1.5
t1[t1>t2]=t2[t1>t2]
M=100

tau=1
x1=pmin(t1, c)
delta1=1*(t1<c)
x2=pmin(t2, c)
delta2=1*(t2<c)

# fit1=dorsurvfun1(x1, delta1, x2, delta2, tau,n.t.grd=1000)
fit=dorfit(x1, delta1, x2, delta2, tau)

fit2 <- dorsurvfun2(x1, delta1, x2, delta2, tau,n.t.grd=500)

plot(c(0, 1), c(0, 1), type="n")

# lines(fit1$time, fit1$surv,col="red")  
# lines(fit$t1, fit$surv1,col="blue") 

lines(fit$t2, fit$surv2) 
lines(fit$t2, fit$surv2m,col="red") 
# lines(fit2$time, fit2$surv,col="green")


#########

