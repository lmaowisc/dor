library(survival)
library(tidyverse)
library(knitr)
## added 
source("DORfunctions.R")

rm(list=ls())
# input: (x1, delta1, x2, delta2, tau)


dorsurvfun1=function(x1, delta1, x2, delta2, tau, n.t.grd=100)
{
  n=length(x1)
  x1=pmin(x1, tau)
  delta1[x1==tau]=1
  x2=pmin(x2, tau)
  delta2[x2==tau]=1
  
  fitc=survfit(Surv(x2, 1-delta2)~1)
  time0=c(0, fitc$time)
  surv0=c(1, fitc$surv)
  

  weight=rep(0, n)
  for(ii in 1:n)
  {if(delta2[ii]>0)
      weight[ii]=1/min(surv0[time0<=x2[ii]])
  }
  
  t.grd=seq(0, tau, length=n.t.grd)
  surv.grd=rep(NA, n.t.grd)
  for(i in 1:n.t.grd)
  { t0=t.grd[i]
    surv.grd[i]=sum((x2>x1+t0)*weight)/n
  }
  
  return(list(time=t.grd, surv=surv.grd))
}


## example

n=200
n.t.grd=100
tau=1.25

B=1000  # number of simulation


surv1 = matrix(0, B, n.t.grd)
surv1LM = matrix(0, B, n.t.grd)
surv1seLM = matrix(0, B, n.t.grd)
surv_cond = matrix(0, B, n.t.grd)
surv_cond_logse = matrix(0, B, n.t.grd)


for(b in 1:B){
  t1=pmin(rexp(n), runif(n, 0, 0.75))
  t2=rexp(n)*1.5
  t1[t1>t2]=t2[t1>t2]
  c=runif(n, 0, 1.5)


  x1=pmin(t1, c)
  delta1=1*(t1<c)
  x2=pmin(t2, c)
  delta2=1*(t2<c)

  fit1=dorsurvfun1(x1, delta1, x2, delta2, tau, n.t.grd=n.t.grd)
  surv1[b,]=fit1$surv
  time0=fit1$time
  ## added by LM
  
  fitLM <- dorfit(x1, delta1, x2, delta2, tau, M=n.t.grd)
  m <- length(fitLM$t1)
  
  ind <- round(ecdf(fitLM$t1)(time0) * m)
  
  surv1LM[b,] <- fitLM$surv1[ind]
  surv1seLM[b,] <- fitLM$se1[ind]
  surv_cond[b,] <- exp(fitLM$log_surv_cond[ind])
  surv_cond_logse[b,] <- fitLM$se_log_surv_cond[ind]
  cat(b,"\n")
}


tmp <- apply(surv_cond, 2, sd)

za <- qnorm(0.975)


results <- tibble(tau = tau, t = time0,
                  surv_LT = colMeans(surv1, na.rm =TRUE),
                  surv_LM = colMeans(surv1LM, na.rm =TRUE),
                  surv_se_LM = colMeans(surv1seLM, na.rm =TRUE),
                  surv_cond = colMeans(surv_cond, na.rm =TRUE),
                  surv_cond_empse = tmp,
                  surv_cond_logse = colMeans(surv_cond_logse, na.rm =TRUE)
                  )




colors <- c("Sepal Width" = "blue", "Petal Length" = "red", "Petal Width" = "orange")

theme_set(theme_bw())

results |> 
  mutate(
    ymin = surv_LM * exp(- za * surv_se_LM/surv_LM),
    ymax = surv_LM * exp(za * surv_se_LM/surv_LM)
  ) |> 
  ggplot(aes(x = t, y = surv_LM)) +
  geom_line(aes(color = "LM"),lwd = 2, alpha = 0.8) +
  geom_line(aes(y = surv_LT, color = "LT"), lwd = 1,  alpha = 0.8) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax),
              alpha = 0.25) +
  coord_cartesian(xlim = c(0, 1.25), ylim = c(0, 1)) +
  scale_color_manual(values = c("LM" = "red", "LT" = "blue")) +
  labs(x = "t",
       y = expression(S[D](t)),
       color = "Legend") 


# Save a single object to a file
saveRDS(results, "results.rds")
# Restore it under a different name
results <- readRDS("results.rds")

tau <- results$tau[1]

par(mfrow=c(1,1))
plot(c(0, tau), c(0, 1), type="n", xlab = "t", ylab = expression(S[D](t)))
lines(time0, apply(surv1, 2, mean), col=2, lwd = 2)  #LT
lines(time0, apply(surv1LM, 2, mean), col=3) #LM
legend("topright", lwd = c(2, 1), col = 2:3, c("LT", "LM"))




n0=100000
t1=pmin(rexp(n0), runif(n0, 0, 0.75))
t2=rexp(n0)*1.5
t1[t1>t2]=t2[t1>t2]


unres <- tibble(dor=c(0, sort((t2-t1)[t2>t1])))

m=length(unres$dor) -1

unres$St <- 1- (0:m)/m

results |> 
  mutate(
    ymin = surv_cond * exp(- za * surv_cond_logse),
    ymax = surv_cond * exp(+ za * surv_cond_logse),
  ) |> 
  ggplot(aes(x = t, y = surv_cond)) +
  geom_line(aes(color = "Restricted"), lwd = 1.5) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3) +
  geom_line(data = unres, aes(x = dor, y = St, color = "Unrestricted"), lwd = 1.5) +
  xlim(0, 1.25) +
  labs(
    x = "Time",
    y = "Survival Probability of DOR",
    color = "Legend"
  ) +
  scale_color_brewer(palette = "Set1")


### inference 
n0=1000000
t1=pmin(rexp(n0), runif(n0, 0, 0.75))
t2=rexp(n0)*1.5
t1[t1>t2]=t2[t1>t2]
t1tau <- pmin(t1, tau)
t2tau <- pmin(t2, tau)
t1tau[t1tau>t2tau]=t2tau[t1tau>t2tau]
res <- tibble(t=c(0, sort((t2tau-t1tau)[t2tau>t1tau])))
m=length(res$t) - 1
res$St <- 1- (0:m)/m

t0 <- c(0.25, 0.5, 0.75, 1)
m <- m + 1

ind <- round(ecdf(res$t)(t0) * m)
res_infer <- res[ind, ]
res_infer$t <- t0


### simulations 
set.seed(2023)
n=200
tau=1.25

B=1000  # number of simulation

m0 <- length(t0)

surv_cond_log = matrix(0, B, m0)
surv_cond_logse = matrix(0, B, m0)
surv_cond_logcp = matrix(0, B, m0)

for(b in 1:B){
  t1=pmin(rexp(n), runif(n, 0, 0.75))
  t2=rexp(n)*1.5
  t1[t1>t2]=t2[t1>t2]
  c=runif(n, 0, 1.5)
  
  
  x1=pmin(t1, c)
  delta1=1*(t1<c)
  x2=pmin(t2, c)
  delta2=1*(t2<c)

  
  fitLM <- dorfit(x1, delta1, x2, delta2, tau, M=n.t.grd)
  m <- length(fitLM$t1)
  
  ind <- round(ecdf(fitLM$t1)(t0) * m)
  
  surv_cond_log[b,] <- fitLM$log_surv_cond[ind]
  surv_cond_logse[b,] <- fitLM$se_log_surv_cond[ind]
  surv_cond_logcp[b,] <- abs(surv_cond_log[b,] - log(res_infer$St)) <= za*surv_cond_logse[b,]
  
  cat(b,"\n")
}

tmp <- apply(exp(surv_cond_log), 2, sd)
results2 <- tibble(t = t0,
                  surv_cond_log = colMeans(surv_cond_log, na.rm =TRUE),
                  surv_cond_empse = tmp,
                  surv_cond_logse = colMeans(surv_cond_logse, na.rm =TRUE),
                  surv_cond_cp = colMeans(surv_cond_logcp, na.rm =TRUE)
)


# Save a single object to a file
saveRDS(results2, "results2.rds")
# Restore it under a different name
results2 <- readRDS("results2.rds")




infer <- merge(res_infer, results2, by = 1)

### tabulate

tab <- infer |> 
  mutate(
    True = St,
    Bias = exp(surv_cond_log) - St,
    SE = surv_cond_empse,
    SEE = exp(surv_cond_log) * surv_cond_logse,
    CP = surv_cond_cp
  ) |> 
  select(t, True, Bias, SE, SEE, CP) 

kable(round(tab, 3), caption = "Estimation and inference of Pr(D > t| D >0).")

