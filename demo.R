library(survival)


# set working directory for DORfunctions.R
setwd("")
source("DORfunctions.R")

## example
############
tau <- 2

n <- 10
##### simulate data ###########
c <- runif(n, 0, 3)
t1 <- rexp(n)
t2 <- rexp(n)*2

t1[t1>t2] <- t2[t1>t2]


x2 <- pmin(t2, c)
delta2 <- (t2<c) + 0 

x1 <- pmin(t1, x2)
delta1 <- (t1< x2) + 0


# fit the estimators of survival functions of D
fit <- dorfit(x1, delta1, x2, delta2, tau)

### naive estimates ####
t1 <- fit$t1
surv1 <- fit$surv1
se1 <- fit$se1

#### improved IPW ####
t2 <- fit$t2
surv2 <- fit$surv2
se2 <- fit$se2
surv2m <- fit$surv2m

# set up the figure
plot(c(0, 2), c(0, 1), type="n",xlab="t",ylab=expression(S[D](t)),cex.lab=1.5,cex.axis=1.5)

# plot the estimated curves
lines(stepfun(t1, c(1,surv1)),col="red",lwd=2) #naive
lines(stepfun(t2, c(1,surv2)),col="blue",lwd=2) #ipw
lines(stepfun(t2, c(1,surv2m)),col="purple",lwd=2) #ipw-monotone

######## get the median ####
# naive
med1 <- t1[which(surv1 <= 0.5)[1]]
# ipw
med2 <- t2[which(surv2 <= 0.5)[1]]
# ipw-mono
med2m <- t2[which(surv2m <= 0.5)[1]]
#########

### area under the curve (mean) ####
# naive
AUC1 <- sum(surv1[-length(surv1)]*diff(t1))
# ipw
AUC2 <- sum(surv2[-length(surv2)]*diff(t2))
# ipw-mono
AUC2m <- sum(surv2m[-length(surv2m)]*diff(t2))


