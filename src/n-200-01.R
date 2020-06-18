# Implement LLSS
# heteroscedastic design
#
# Date          Changed
#----------------------------------------------------------------------------
# 2016-07-22     Original Code
# 2016-08-26     objective function of tau.hat corrected at each step.

rm(list=ls())
library("xtable")
library("lattice")
library(quantreg)
start.time = proc.time()[3]

source('C-functions.R')

# PARAMETERS
m.quant=0.5
seed.n="01"
n=200
set.seed(as.numeric(seed.n))
ap.0 = ap.a + (qnorm(m.quant,mean=0,sd=std)) * het.0
T.true= support(ap.0,tr)
tau.0 = 0.5
const= rep(1, n)

nrep=5
alpha= 0.1
c=1.1
signal= 1
std= 0.5
p=500
tr= 10^(-6)
p.2=p/2
bt = c(0,qnorm(0.75,mean=0,sd=std),0,0,rep(0,p.2-4))
#bt = c(0,0,0,0,rep(0,p.2-4))
dt = c(0,1,0,0,rep(0,p.2-4))
ap.a = c(bt,dt)
het.0 = c(0,1,0,rep(0,p.2-3),0,0,0,rep(0,p.2-3))

f.grid<-function(q,low,high){
  sort.q=sort(q)
  n.q=length(q)
  sort.q[(low*n.q):(high*n.q)]
}

Sigma= toeplitz((1/2)^{0:(p.2-2)})
Sigma.sqroot= chol(Sigma)



cat('model quantile=',m.quant,'\n')
cat('c=',c,'\n')
cat('alpha=',alpha,'\n')
cat('std=',std,'\n')
cat('n=',n,'\n')
cat('p=',p,'\n')
cat('nrep=',nrep,'\n')

# MONTE-CARLO
ap.hat.o1 <- matrix(0, nrep, p)  # ap.hat of oracle 1
ap.hat.o2 <- matrix(0, nrep, p)  # ap.hat of oracle 2
tau.hat.o2 = matrix(0,nrep,1) # tau.hat of oracle 2

ap.hat.1 = matrix(0,nrep,p) # alpha-hat for step 1
tau.hat.1 = matrix(0,nrep,1) # tau-hat for step 1
tau.hat.2 = matrix(0,nrep,1) # tau-hat for step 2
tau.hat.3a = matrix(0,nrep,1) # tau-hat for step 3a
tau.hat.3b = matrix(0,nrep,1) # tau-hat for step 3a
ap.hat.3a = matrix(0,nrep,p)
ap.hat.3b = matrix(0,nrep,p)
M= matrix(0, nrep, 1)
K= matrix(0, nrep, 1)
M3= matrix(0, nrep, 1)
K3= matrix(0, nrep, 1)
M3b= matrix(0, nrep, 1)
K3b= matrix(0, nrep, 1)

CI.o2 = matrix(0,nrep,2)
CI.1 = matrix(0,nrep,2) # confidence intervals
CI.2 = matrix(0,nrep,2) # confidence intervals
CI.3a = matrix(0,nrep,2) # confidence intervals
CI.3b = matrix(0,nrep,2) # confidence intervals



for(j in 1:nrep){
  cat('----------------------------','\n')
  cat('replication: ',j,'\n')
  #--------------------------------------------------------------
  # DGP
  #--------------------------------------------------------------
  dat <- dgp(n,p,p.2,const,het.0,ap.a,tau.0)
  Z= dat$Z
  X.2= dat$X.2
  Q = dat$Q
  Y= dat$Y
  grid.tau = f.grid(Q,0.15,0.85)
  grid.size = length(grid.tau)

  #-----------m.quant---------------------------------------------------
  # Oracle 1: J_0 & tau_0 known. Estimate ap
  #--------------------------------------------------------------
  T.oracle<- T.true
  X.oracle<- cbind(X.2,X.2*(Q>tau.0))[,T.oracle]
  ap.hat.o1[j,T.oracle]<- rq(Y~-1+X.oracle, tau=m.quant)$coef


  #--------------------------------------------------------------
  # Oracle 2: J_0 known. Estimate ap & tau
  #--------------------------------------------------------------
  obj.o2 = rep(0,grid.size)
  ap.0.tilde.grid = matrix(0,grid.size,p)
  for (k in (1:grid.size)) {
    ind = (Q > grid.tau[k])
    X.sel = cbind(const,Z)
    X.sel = cbind(X.sel, X.sel*ind)
    X.sel<- X.sel[,T.true]
    m.o2 = rq(Y~-1+X.sel, tau=m.quant)
    ap.0.tilde.grid[k,T.true] = m.o2$coef
    u.hat = m.o2$resid
    obj.o2[k]= sum(check.fn(u=u.hat,tau=m.quant))
  }
  pick.min = which.min(obj.o2)
  if (pick.min==1 | pick.min==grid.size ) {
    tau.hat.o2[j] = grid.tau[pick.min]} # boundary solution. pick as it is
  else {
    tau.hat.o2[j] = (grid.tau[pick.min]+grid.tau[pick.min+1])/2 # otherwise, pick the mid-point
  }
  ap.hat.o2[j,] <- ap.0.tilde.grid[pick.min,]

  CI.o2[j,]=tau.inference(ap.hat.o2[j,],tau.hat.o2[j,1],bar.T=n*0.5,X.for.bt=X.2,X.for.dt=X.2,dt.hat=ap.hat.o2[j,251:500],m.quant=m.quant)



  #--------------------------------------------------------------
  # Step 1
  #--------------------------------------------------------------
  ap.hat.grid = matrix(0,grid.size,p)
  obj.s1 = matrix(0,grid.size,1)
  # choose rq.lambda as a sup over tau.grid
  rq.lambda.grid = matrix(0,grid.size,1)
  for (i.tau in 1:grid.size){
    X = cbind(X.2,X.2*(Q>grid.tau[i.tau])) # inequality change
    rq.lambda.grid[i.tau,] = Lambda(X[,-1], R.rep = 1000, tau=m.quant , c=c, alpha)
  }
  max.rq.lambda = max(rq.lambda.grid)
  for (i.tau in 1:grid.size){
    X = cbind(X.2,X.2*(Q>grid.tau[i.tau])) # inequality change
    rq.lambda = max.rq.lambda*c(1,apply(X[,-1],2,norm2n))
    fit= rq(Y ~ X[,-1], tau=m.quant, method="lasso",lambda = rq.lambda)
    ap.hat.grid[i.tau,] = fit$coef
    obj.s1[i.tau,1] = fit$rho + sum(rq.lambda * abs(ap.hat.grid[i.tau,]))
  }
  opt = which.min(obj.s1)
  if (opt==1 | opt==grid.size ) {
    tau.hat.1[j,1]=grid.tau[opt] # boundary solution. pick as it is
  }
  else {
    tau.hat.1[j,1] = (grid.tau[opt]+grid.tau[opt+1])/2 # otherwise, pick the mid-point
  }
  ap.hat.1[j,]=ap.hat.grid[opt,]

  T.hat = support(ap.hat.1[j,],tr)
  m.hat = length(setdiff(T.hat, T.true))
  k.hat= length(setdiff(T.true, T.hat))

  M[j,]= m.hat
  K[j,]= k.hat

  CI.1[j,]=tau.inference(ap.hat.1[j,],tau.hat.1[j,1],bar.T=n*0.5,X.for.bt=X.2,X.for.dt=X.2,dt.hat=ap.hat.1[j,251:500],m.quant=m.quant)


  #--------------------------------------------------------------
  # Step 2
  #--------------------------------------------------------------
  obj.s2= matrix(0,grid.size,1)
  for (i.tau in 1:grid.size){
    X = cbind(X.2,X.2*(Q>grid.tau[i.tau])) # inequality change
    Y.hat = X %*% ap.hat.1[j,]
    obj.s2[i.tau,1] = sum(check.fn(Y-Y.hat, tau=m.quant))
  }
  opt.2 = which.min(obj.s2)
  if (opt.2==1 | opt.2==grid.size ) {
    tau.hat.2[j,1]=grid.tau[opt.2] # boundary solution. pick as it is
  }
  else {
    tau.hat.2[j,1] = (grid.tau[opt.2]+grid.tau[opt.2+1])/2 # otherwise, pick the mid-point
  }

  CI.2[j,]=tau.inference(ap.hat.1[j,],tau.hat.2[j,1],bar.T=n*0.5,X.for.bt=X.2,X.for.dt=X.2,dt.hat=ap.hat.1[j,251:500],m.quant=m.quant)

  #--------------------------------------------------------------
  # Step 3a
  #--------------------------------------------------------------
  X.sel = cbind(X.2,X.2*(Q>tau.hat.2[j])) # inequality change
  rq.lambda.3 = lambda.BC(X.sel[,-1], tau=m.quant, c=c,alpha=alpha)
  fit.3= rq(Y ~ X.sel[,-1], tau=m.quant ,method='lasso',lambda=rq.lambda.3)
  ap.hat.3a[j,] = fit.3$coef

  T.hat3 = support(ap.hat.3a[j,],tr)
  m.hat = length(setdiff(T.hat3, T.true))
  k.hat= length(setdiff(T.true, T.hat3))

  M3[j,]= m.hat
  K3[j,]= k.hat

  # Update tau.hat
  obj.3a= matrix(0,grid.size,1)
  for (i.tau in 1:grid.size){
    X = cbind(X.2,X.2*(Q>grid.tau[i.tau])) # inequality change
    Y.hat = X %*% ap.hat.3a[j,]
    obj.3a[i.tau,1] = sum(check.fn(Y-Y.hat, tau=m.quant))
  }
  opt.3a = which.min(obj.3a)
  if (opt.3a==1 | opt.3a==grid.size ) {
    tau.hat.3a[j,1]=grid.tau[opt.3a] # boundary solution. pick as it is
  }
  else {
    tau.hat.3a[j,1]= (grid.tau[opt.3a]+grid.tau[opt.3a+1])/2 # otherwise, pick the mid-point
  }

  CI.3a[j,]=tau.inference(ap.hat.3a[j,],tau.hat.3a[j,1],bar.T=n*0.5,X.for.bt=X.2,X.for.dt=X.2,dt.hat=ap.hat.3a[j,251:500],m.quant=m.quant)

  #--------------------------------------------------------------
  # Step 3b
  #--------------------------------------------------------------
  X.sel.b = cbind(X.2,X.2*(Q>tau.hat.2[j])) # inequality change
  mu.n = rq.lambda.3[1]/n
  w = (abs(ap.hat.3a[j,]) < mu.n) + (abs(ap.hat.3a[j,]) >= mu.n)*(abs(ap.hat.3a[j,]) < 3.7*mu.n)*((3.7*mu.n-abs(ap.hat.3a[j,]))/(mu.n*2.7))
  sigs = apply(X.sel.b,2,norm2n)
  rq.lambda.3b = as.numeric((log(log(n))*n*mu.n)*w*sigs)
  fit.3b= rq(Y ~ X.sel.b[,-1], tau=m.quant,method='lasso',lambda=rq.lambda.3b)
  ap.hat.3b[j,] = fit.3b$coef

  T.hat3b = support(ap.hat.3b[j,],tr)
  m.hat = length(setdiff(T.hat3b, T.true))
  k.hat= length(setdiff(T.true, T.hat3b))

  M3b[j,]= m.hat
  K3b[j,]= k.hat

  # Update tau.hat
  obj.3b= matrix(0,grid.size,1)
  for (i.tau in 1:grid.size){
    X = cbind(X.2,X.2*(Q>grid.tau[i.tau])) # inequality change
    Y.hat = X %*% ap.hat.3b[j,]
    obj.3b[i.tau,1] = sum(check.fn(Y-Y.hat, tau=m.quant))
  }
  opt.3b = which.min(obj.3b)
  if (opt.3a==1 | opt.3a==grid.size ) {
    tau.hat.3b[j,1]=grid.tau[opt.3b] # boundary solution. pick as it is
  }
  else {
    tau.hat.3b[j,1]= (grid.tau[opt.3b]+grid.tau[opt.3b+1])/2 # otherwise, pick the mid-point
  }

  CI.3b[j,]=tau.inference(ap.hat.3b[j,],tau.hat.3b[j,1],bar.T=n*0.5,X.for.bt=X.2,X.for.dt=X.2,dt.hat=ap.hat.3b[j,251:500],m.quant=m.quant)

  print.mat = matrix(NA,6,length(ap.0[T.true])+1)
  rownames(print.mat)=c('ap.0','ap.o1','ap.o2','ap.1','ap.3a','ap.3b')
  colnames(print.mat)=c(T.true,'tau')

  print.mat[1,] = c(ap.0[T.true],tau.0)
  print.mat[2,] = c(ap.hat.o1[j,T.true],tau.0)
  print.mat[3,] = c(ap.hat.o2[j,T.true],tau.hat.o2[j])
  print.mat[4,] = c(ap.hat.1[j,T.true],tau.hat.1[j])
  print.mat[5,] = c(ap.hat.3a[j,T.true],tau.hat.3a[j])
  print.mat[6,] = c(ap.hat.3b[j,T.true],tau.hat.3b[j])
  print(round(print.mat,4))


  cat('ap.1.non.zero=',support(ap.hat.1[j,],tr),'\n')
  cat('ap.3a.non.zero=',support(ap.hat.3a[j,],tr),'\n')
  cat('ap.3b.non.zero=',support(ap.hat.3b[j,],tr),'\n')


}

ap.hat.1=ap.hat.1*(abs(ap.hat.1)>tr)
ap.hat.3a=ap.hat.3a*(abs(ap.hat.3a)>tr)
ap.hat.3b=ap.hat.3b*(abs(ap.hat.3b)>tr)


# Print computation time & save the result
runt = proc.time()[3]-start.time
runt_h = floor(runt/3600)
runt_m = floor( (runt/3600 - runt_h) * 60 )
runt_s = floor( (runt/60 - (runt_h*60) - runt_m) * 60 )
cat('\n',"runtime = ",runt_h,"h",runt_m,"m",runt_s,"s",'\n','\n')
cat('-----------------------------------------------------------','\n')

save.image(file =paste("../results/n-",n,'-',seed.n,".RData",sep=''))
