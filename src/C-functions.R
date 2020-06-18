# PARAMETERS
nrep=1
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


# FUNCTIONS

dgp <- function(n,p,p.2,const,het.0,ap.a,tau.0){
  Z= cbind(matrix(rnorm(n*(p.2-1), mean=5, sd=1), n, (p.2-1)  ))%*%Sigma.sqroot
  #Z= cbind(matrix(rchisq(n*(p.2-1), df=3), n, (p.2-1)  ))/(2*sqrt(6))
  X.2= cbind(const, Z )
  Q = runif(n)
  #U = (runif(n)-0.5)
  U = rnorm(n, mean=0, sd=std)
  u.a =  het.0 %*% t(U)
  ap = matrix(ap.a,p,n) + u.a

  Y= diag(cbind(X.2, X.2*(Q>tau.0)) %*% ap)
  return(list(Z=Z,X.2=X.2,Q=Q,U=U,Y=Y))
}




norm2n= function(z){  sqrt(mean((z)^2)) }
#norm2n= function(z){  sqrt(mean((z-mean(z))^2)) }

norm2 = function(z){sqrt(sum(z^2))}
lambda.BC= function(X, R.rep = 1000, tau, c, alpha){
  n = nrow(X)
  sigs = apply(X,2,norm2n)
  U = matrix(runif(n * R.rep),n)
  R.mat = (t(X) %*% (tau - (U < tau)))/(sigs*sqrt(tau*(1-tau)))  # 500 x 1000 division row-wise
  r = apply(abs(R.mat),2,max)
  c * quantile(r, 1 - alpha) * sqrt(tau*(1-tau))*c(1,sigs)
}

Lambda= function(X, R.rep = 1000, tau, c, alpha){
  n = nrow(X)
  sigs = apply(X,2,norm2n)
  U = matrix(runif(n * R.rep),n)
  #R.mat = (t(X) %*% (tau - (U < tau)))/(sigs)  # 500 x 1000 division row-wise
  R.mat = (t(X) %*% (tau - (U < tau)))/(sigs*sqrt(tau*(1-tau)))
  r = apply(abs(R.mat),2,max) # 1000x1 simulated values
  #c * quantile(r, 1 - alpha)
  c * quantile(r, 1 - alpha) * sqrt(tau*(1-tau))
}


support= function(x, tr) {
  m= rep(0, length(x))
  for (i in 1:length(x)) if( abs(x[i])> tr ) m[i]= i
  m = m[m>0]
  m
}

er = function(bt.hat,Sig,bt.0){
  rep = nrow(bt.hat)
  n.par = length(bt.0)
  G = bt.hat-matrix(bt.0,rep,n.par,byrow=T)
  M = G %*% Sig %*% t(G)
  N = diag(M)
  O = mean (N)
  sqrt(O)
}


er.sim = function(bt.hat,bt.0,tau.hat,tau.0,sim.size){
  rep.size = length(tau.hat)
  er.sim.dat <- dgp(sim.size,p,p.2,const,het.0,ap.a,tau.0)
  q = er.sim.dat$Q
  z = er.sim.dat$Z
  x.2 = cbind(1,z)
  x = cbind(x.2, x.2*(q>tau.0))
  ind.0 = x %*% bt.0

  ms = matrix(0,rep.size,1)
  for (i in (1:rep.size)){
    x.hat = cbind(x.2, x.2*(q>tau.hat[i])) # inequality change
    ind.hat = x.hat %*% bt.hat[i,]
    ms[i] = mean((ind.hat-ind.0)^2)
  }
  mean(sqrt(ms)) # switch mean and sqrt. mean prediction error over replications
}


l2norm.bias = function(bt.hat,bt.0){
  means= (apply(bt.hat, 2, mean)-bt.0)
  sqrt(sum(means^2))
}

loss.fun <- function(X,Y,Q,tau,ap,q){
  X.tau=cbind(X, X*(Q>tau)) # inequality change
  index=X.tau %*% ap
  U=Y-index
  U*(q-(U <= 0))
}


poisson.p<-function(ld,max.T){
  seq.t = 0
  t = 0
  repeat{
    U=runif(1)
    t=t+(-(1/ld)*log(U))
    if (t>max.T) break
    seq.t=c(seq.t,t)
  }
  return(seq.t)
}

# Count the number of events before [tt] given [seq.t]
count.N <-function(seq.t,tt){
  sum(seq.t<tt)-1
}

check.fn=function(u,tau){ u*(tau-(u<0)) }


excess.risk.sim = function(ap.hat,ap.0,tau.hat,tau.0,sim.size,qt){
  rep.size = length(tau.hat)
  excess.risk=rep(0,rep.size)

  excess.risk.sim.dat <- dgp(sim.size,p,p.2,const,het.0,ap.a,tau.0)
  q = excess.risk.sim.dat$Q
  z = excess.risk.sim.dat$Z
  x.2 = cbind(1,z)
  y= excess.risk.sim.dat$Y

  for (i in (1:rep.size)) {
    # Calculate two loss functions with estimated parameters and true parameters
    loss.hat <- loss.fun(x.2,y,q,tau.hat[i],ap.hat[i,],q=qt)
    loss.true <- loss.fun(x.2,y,q,tau.0,ap.0,q=qt)
    excess.risk[i] <-  mean(loss.hat-loss.true)
  }
  mean(excess.risk)
}


tau.inference <- function(ap.hat, tau.hat, bar.T, X.for.bt, X.for.dt, dt.hat, m.quant){
  B=1000
  #B=1
  h.B = rep(NA,B)

  for (k in (1:B)){
    # Estimate the jump rate f_Q(tau_0)
    f.Q=approxfun(density(Q,bw='nrd')) # bw choice: 1.06*min(s,interquartile-range/1.34)*n^(-1/5) rule-of-thumb
    jump.rate=f.Q(tau.hat)
    jump.h1 = poisson.p(ld=jump.rate,max.T=bar.T)
    jump.h2 = poisson.p(ld=jump.rate,max.T=bar.T)

    max.N1 = length(jump.h1[-1])
    max.N2 = length(jump.h2[-1])
    # Simulate rho.1i and rho.2i
    X.all = as.matrix(cbind(X.for.bt,X.for.dt*(Q>tau.hat)))
    Y.hat = X.all %*% ap.hat
    eps.i = Y-Y.hat
    e.rho.1i=check.fn(eps.i-X.for.dt%*%dt.hat,tau=m.quant)-check.fn(eps.i,tau=m.quant)
    e.rho.2i=check.fn(eps.i+X.for.dt%*%dt.hat,tau=m.quant)-check.fn(eps.i,tau=m.quant)

    cdf.rho1=ecdf(e.rho.1i)
    cdf.rho2=ecdf(e.rho.2i)
    rho.1i=quantile(cdf.rho1,runif(max.N1))
    rho.2i=quantile(cdf.rho2,runif(max.N2))

    # Grid search for h1 and h2, respectively
    M.h.grid = rep(0,max.N1+max.N2)
    h.grid = c(-jump.h1[-1],jump.h2[-1])
    M.h.grid[1] = rho.1i[1]
    for (i in (2:max.N1)){
      M.h.grid[i]=M.h.grid[i-1] + rho.1i[i]
    }
    M.h.grid[max.N1+1]= rho.2i[1]
    for (i in ((max.N1+2):(max.N1+max.N2))){
      M.h.grid[i]=M.h.grid[i-1] + rho.2i[i-max.N1]
    }

    opt.M.h = which.min(M.h.grid)
    #plot(h.grid,M.h.grid)
    opt.h = h.grid[opt.M.h]
    h.B[k]=opt.h
  }

  h.B.n = h.B/n
  l.bnd = quantile(h.B.n,0.025)
  u.bnd = quantile(h.B.n,0.975)
  cat('l.bnd=',l.bnd,'\n')
  cat('u.bnd=',u.bnd,'\n')
  CI.tau = c(tau.hat+l.bnd, tau.hat+u.bnd)
  cat('95% confidence interval for tau','\n')
  cat('[',CI.tau,']','\n')
  return(CI.tau)
}


oracle.prop<-function(true.index,ap.hat.s){
  total = nrow(ap.hat.s)
  oracle.mat = ap.hat.s
  oracle.mat[,true.index]=(ap.hat.s[,true.index]!=0)
  oracle.mat[,-true.index]=(ap.hat.s[,-true.index]==0)
  true.detect=sum(apply(oracle.mat,1,prod))
  o.p = true.detect / total
  return(o.p)
}

mse <- function(true.matrix, ap.hat.s){
  sq.er = (ap.hat.s-true.matrix)^2
  return(apply(sq.er,2,mean))
}
