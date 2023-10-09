library(survival)

rm(list=ls())
# input: (x1, delta1, x2, delta2, tau)

dorsurvfun2=function(x1, delta1, x2, delta2, tau, n.t.grd=100)
{
n=length(x1)
x1=pmin(x1, tau)
delta1[x1==tau]=1
x2=pmin(x2, tau)
delta2[x2==tau]=1

fitc=survfit(Surv(x2, 1-delta2)~1)
time0=c(0, fitc$time)
surv0=c(1, fitc$surv)


t.grd=seq(0, tau, length=n.t.grd)
surv.grd=rep(NA, n.t.grd)
for(i in 1:n.t.grd)
{ t0=t.grd[i]
delta=(delta2+(1-delta2)*delta1)*(x2>x1+t0)     
tt=pmin(x2, x1+t0, tau)
prob=rep(NA, n)
for(ii in 1:n)
{prob[ii]=min(surv0[time0<=tt[ii]])}
surv.grd[i]=mean(delta/prob)
}

return(list(time=t.grd, surv=surv.grd))
}



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
  
  tt=x2
  prob=rep(NA, n)
  for(ii in 1:n)
      {prob[ii]=min(surv0[time0<=tt[ii]])}
  
  t.grd=seq(0, tau, length=n.t.grd)
  surv.grd=rep(NA, n.t.grd)
  for(i in 1:n.t.grd)
  { t0=t.grd[i]
    delta=delta2*(x2>x1+t0)     
    surv.grd[i]=mean(delta/prob)
  }
  
  return(list(time=t.grd, surv=surv.grd))
}


## example

n=100
c=runif(n, 0, 1.5)
t1=rexp(n)
t2=rexp(n)*1.5
t1[t1>t2]=t2
n.t.grd=100

tau=1
x1=pmin(t1, c)
delta1=1*(t1<c)
x2=pmin(t2, c)
delta2=1*(t2<c)

fit1=dorsurvfun1(x1, delta1, x2, delta2, tau, n.t.grd=100)
fit2=dorsurvfun2(x1, delta1, x2, delta2, tau, n.t.grd=100)

plot(c(0, 1), c(0, 1), type="n")
lines(fit1$time, fit1$surv)  
lines(fit2$time, fit2$surv)


