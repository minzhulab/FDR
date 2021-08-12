
## plot mixture normal component
plot_mix_comps <- function(x, mu, sigma, pi) {   
  pi * dnorm(x, mu, sigma)
}


### rejection region for given lambda
# Significance region Re(lambda):  Story (2003) 
Re <- function(lambda, mu0=0, s0=1, mu1,s1, pi0)
{
  c1 = (mu0/s0)^2 - (mu1/s1)^2
  c2 = 2*log((1-lambda)/lambda*pi0/(1-pi0)*s1/s0)
  
  a = 1/s0^2 - 1/s1^2
  b = -2*mu0/s0^2 + 2*mu1/s1^2
  c = c1-c2 
  #browser() 
  sol = NA
  if(a == 0)  sol = -c/b
  if(a !=0 )  sol = c( (-b-sqrt(b^2-4*a*c))/(2*a), (-b+sqrt(b^2-4*a*c))/(2*a))
  return(list(sol= sort(sol), a=a))
}


## theoretical FDR FNR calculation
FDRFNR <- function(par, mu0 = 0, s0 = 1, mu1, s1, pi0)
{
  if(par$a == 0)  {
    th = par$sol[1]
    k1 = pi0*pnorm(th,  mean=mu0, sd=s0, lower.tail = FALSE)  # pi0*Pr(T in R|H=0)
    k2 = (1-pi0)*pnorm(th,mean=mu1, sd=s1)   # pi1*Pr(T not in R|H=1)
    
    FDR = k1/(k1+(1-pi0)*pnorm(th,mean=mu1, sd=s1, lower.tail = FALSE))
    FNR = k2/(pi0*pnorm(th,  mean=mu0, sd=s0)+k2)
    #be= w*FDR/pi0 + (1-w)*FNR/(1-pi0)
  }
  else if(par$a > 0) {
    th = par$sol
    k1 = pi0*(pnorm(th[1],  mean=mu0, sd=s0) + pnorm(th[2],  mean=mu0, sd=s0, lower.tail = FALSE))  # pi0*Pr(T in R|H=0)
    k2 = (1-pi0)*(pnorm(th[2],  mean=mu1, sd=s1) - pnorm(th[1],  mean=mu1, sd=s1))   # pi1*Pr(T not in R|H=1)
    k3 = (1-pi0)*(pnorm(th[1],  mean=mu1, sd=s1) + pnorm(th[2],  mean=mu1, sd=s1, lower.tail = FALSE))   # pi1*Pr(T in R | H=1)
    k4 = pi0*(pnorm(th[2],  mean=mu0, sd=s0) - pnorm(th[1],  mean=mu0, sd=s0))  # pi0*Pr(T not in R | H=0)
    FDR =  k1/(k1+k3)    
    FNR = k2/(k4+k2)

  }
  else {
    th =  par$sol
    k1 = pi0*(pnorm(th[2],  mean=mu0, sd=s0) - pnorm(th[1],  mean=mu0, sd=s0))   # pi0*Pr(T in R|H=0)
    k2 = (1-pi0)*(pnorm(th[1],  mean=mu1, sd=s1)+ pnorm(th[2],  mean=mu1, sd=s1, lower.tail = FALSE))   # pi1*Pr(T not in R|H=1)
    k3 = (1-pi0)*(pnorm(th[2],  mean=mu1, sd=s1)- pnorm(th[1],  mean=mu1, sd=s1))     # pi1*Pr(T in R | H=1)
    k4 = pi0*(pnorm(th[1],  mean=mu0, sd=s0) + pnorm(th[2],  mean=mu1, sd=s1, lower.tail = FALSE))   # pi0*Pr(T not in R | H=0)
    FDR =  k1/(k1+k3)     
    FNR = k2/(k4+k2)
  }
  return(list(FDR=FDR, FNR=FNR))
}



# Bayes error/loss 
BE <- function(lambda, w, mu0=0, s0=1, mu1, s1, pi0, t.df=NA, lossonly=TRUE)
{
  th = Re(lambda, mu0= mu0, s0 = s0, mu1= mu1, s1=s1, pi0=pi0)
  error = FDRFNR(par=th, mu0= mu0, s0 = s0, mu1= mu1, s1=s1,pi0=pi0)
  be= w*error$FDR + (1-w)*error$FNR
  
  if(lossonly) return (be)   # Bayes error
  else {
    cz = min(th$sol[th$sol>0])
    if(is.na(t.df)) thurdle = qnorm(0.5+pnorm(cz)/2)  else thurdle = qt(0.5+pnorm(cz)/2, df= t.df)
    return(list(loss=be, FDR = error$FDR, FNR = error$FNR, thurdle = thurdle, rej.region =th$sol, a=th$a))
  }
}


## target a pre-given FDR
targetFDR <- function(FDR=0.05,mu0=0, s0=1, mu1,s1, pi0)
{
  fnn <- function(lambda,mu0, s0, mu1,s1, pi0)
  {
    th = Re(lambda, mu0= mu0, s0 = s0, mu1= mu1, s1=s1, pi0=pi0)
    error = FDRFNR(par=th, mu0= mu0, s0 = s0, mu1= mu1, s1=s1,pi0=pi0)
    return( (error$FDR - FDR)^2 )
  }

  obj = optim(par = 0.1, mu0= mu0, s0 = s0, mu1= mu1, s1=s1, pi0=pi0, fn = fnn, method="L-BFGS-B", lower=0, upper=0.8, control = list(maxit = 1000))  
  #browser()
  lambda1 = obj$par
  lambda1
}

### emprical FDR FNR calculation based on posterior prob
FDRTable <- function(tx, xseq)
{
  colnames <- c("N00","N01","N10","N11","fdr","fpr","fnegr","fnr","sens","spec")
  xm <- matrix(nrow=length(xseq),ncol=10,
               dimnames=list(paste('c0=',xseq,sep=''), colnames))
  xmrow <- 1
  for (c0 in xseq)
  {
    a <- sum (tx[tx > c0])
    b <- sum (tx[tx <= c0])
    c <- sum ( 1-tx[tx > c0])
    d <- sum ( 1-tx[tx <= c0])
    fpr <- b/(a+b)
    fnr <- c/(c+d)
    
    fdr <- b/(b+d)  
    
    fndr <- c/(a+c)
    sens <- d/(c+d)
    spec <- a/(a+b)
    
    xm[xmrow,] <- c(a,b,c,d,fdr,fpr,fnr,fndr,sens,spec)
    xmrow <- xmrow + 1
  }
  xm
}


## for p-hacking
twomomentmatch <- function(p0, a, b,z)  # z is two-component normal mixture,  a, b are truncation points 
{
  obs = c(mean(z), mean(z^2))   # observed first two sample moments
  # sn + normal mixture
  ini = c(2, log(2))
  obj = optimx(par=ini, fn = mom.match, control = list(maxit = 1000), a=a, b=b, pi0=p0, target=obs)    #   "BFGS", "Nelder-Mead"
  
  m1 = obj$p1[which.min(obj$value)]
  s1 = exp(obj$p2[which.min(obj$value)])
  
  return(c(m1, s1))
}

phackingtab <- function(est, x, FDR=0.05){
  res = data.frame(est)
  res$thurdle = 0
  res$FNR = 0

  for(i in 1:nrow(est))
  {
    lamb <- targetFDR(FDR=FDR,mu1=est[i,2], s1=est[i,3], pi0=est[i,1])   # lambda which is associated with a given FDR level
    output = BE(lambda=lamb, w=1/2, mu1=est[i,2], s1=est[i,3], pi0=est[i,1], lossonly=F)
    res$thurdle[i] = output$thurdle
    res$FNR[i] = output$FNR
  }
  names(res)[1:3] <- c("pi0", "mu1","sigma1")
  res
} 


###  first two moments for left and two-sided truncation
tr.moment <- function(mu, sigma, a, b=Inf)
{
  if(is.infinite(b)){  #left-sided truncation
    theta = (a-mu)/sigma
    mom1 = mu + sigma*dnorm(theta)/(1-pnorm(theta))
    mom2 = mu^2+sigma^2+(2*mu*sigma + theta*sigma^2)*dnorm(theta)/(1-pnorm(theta))
  } else{  # two-sided truncation
    theta1 = (a-mu)/sigma
    theta2 = (b-mu)/sigma
    cc1= (dnorm(theta1)-dnorm(theta2))/(pnorm(theta2)-pnorm(theta1))
    cc2 = (theta1*dnorm(theta1) - theta2*dnorm(theta2))/(pnorm(theta2)-pnorm(theta1))
    mom1 = mu + sigma*cc1
    mom2 = mu^2+sigma^2+sigma^2*cc2 + 2*mu*sigma*cc1
  }
  
  return(c(mom1,mom2))
}


## first two moments for a two-component normal dist (truncated)  
mom.match <- function(para, a, b=Inf, pi0, target) 
{
  if(length(para) == 2){
    m0 = 0
    s0 = 1   
    m1 = para[1]
    s1 = exp(para[2])  
  } 
  if(length(para) == 4){
    m0 = para[1]
    s0 = exp(para[2])   
    m1 = para[3]
    s1 = exp(para[4])  
  }
  
  P0 = pnorm(b, mean=m0, sd=s0) - pnorm(a, mean=m0, sd=s0)
  P1 = pnorm(b, mean=m1, sd=s1) - pnorm(a, mean=m1, sd=s1)
  newpi0 = pi0*P0/(pi0*P0 + (1-pi0)*P1)
  
  com1 = tr.moment(mu=m0, sigma=s0, a=a, b=b)
  com2 = tr.moment(mu=m1, sigma=s1, a=a, b=b)
  com = newpi0*com1 + (1-newpi0)*com2
  
  sum((com-target)^2) 
}


# ###### loglikelihood for truncated normal
tr.loglik.mixnormal <- function(pi0, mu0=0, sigma0=1, mu1, sigma1, x, a, b=Inf)
{
  P0 = pnorm(b, mean=mu0, sd=sigma0) - pnorm(a, mean=mu0, sd=sigma0)
  P1 = pnorm(b, mean=mu1, sd=sigma1) - pnorm(a, mean=mu1, sd=sigma1)
  newpi0 = pi0*P0/(pi0*P0 + (1-pi0)*P1)

  ans = newpi0*dnorm(x, mean=mu0, sd=sigma0)/P0 + (1-newpi0)*dnorm(x, mean=mu1, sd=sigma1)/P1
  return(list(newpi0= newpi0, loglik = sum(log(ans))))
}



