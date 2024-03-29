library(survival)

source("DORfunctions.R")

## example

n=200
# set restriction time
tau=1.25
# simulate data
set.seed(2023)
t1=pmin(rexp(n), runif(n, 0, 0.75))
t2=rexp(n)*1.5
t1[t1>t2]=t2[t1>t2]
c=runif(n, 0, 1.5)


x1=pmin(t1, c)
delta1=1*(t1<c)
x2=pmin(t2, c)
delta2=1*(t2<c)


# fit the data
fitLM <- dorfit(x1, delta1, x2, delta2, tau, med_inf = TRUE)

za <- qnorm(0.975)
## marginal survival function P(D > t)
t <- fitLM$t1 # times
St <- fitLM$surv1 # P(D > t)
se <- fitLM$se1 # se of P(D > t)

## conditional survival function P(D > t|D > 0)
cSt <- exp(fitLM$log_surv_cond)
cSt_U <- exp(fitLM$log_surv_cond + za * fitLM$se_log_surv_cond) # 
cSt_L <- exp(fitLM$log_surv_cond - za * fitLM$se_log_surv_cond)


par(mfrow=c(1, 1))

# set up the figure
plot(c(0, 1.2), c(0, 1), type="n",xlab="t",ylab=expression(S[D](t)),cex.lab=1.2,cex.axis=1.2)

# plot the marginal curve P(D > t)
lines(stepfun(t, c(1,St)),col="red",lwd=2, do.points = FALSE) # marginal
# lines(stepfun(t, c(1,St + za * se)),col="red",lwd=1, lty = 3, do.points = FALSE) # upper 95% CI 
# lines(stepfun(t, c(1,St - za * se)),col="red",lwd=1, lty = 3, do.points = FALSE) # lower 95% CI 

# plot the conditional curve P(D > t|D > 0)
lines(stepfun(t, c(1,cSt)),col="blue",lwd=2, do.points = FALSE) # conditional
# lines(stepfun(t, c(1,cSt_U)),col="blue",lwd=1, lty = 3, do.points = FALSE) # upper 95% CI 
# lines(stepfun(t, c(1,cSt_L)),col="blue",lwd=1, lty = 3, do.points = FALSE) # lower 95% CI 

legend("topright", lwd = 2, col = c("red", "blue"), c("P(D>t)", "P(D>t|D>0)"))

######## get the median ####
med <- fitLM$med
med
# [1] 0.551916

med_lo <- fitLM$med_lo
med_hi <- fitLM$med_hi
med_lo
# [1] 0.4018724
med_hi
# [1] 0.7019595

##### get the AUC

rmst <- fitLM$rmst
rmst

rmst_lo <- fitLM$rmst_lo
rmst_hi <- fitLM$rmst_hi
rmst_lo
# [1,] 0.4862219
rmst_hi
# [1,] 0.664599
