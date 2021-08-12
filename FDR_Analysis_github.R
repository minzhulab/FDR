
#### here is how to install qvalue package
# install.packages("BiocManager")
# BiocManager::install("qvalue")


load("datasets.RData")   #load data sets
ls()

################ Data sets description  ################
## dataset "factortv" contains t-values for published anomaly factors from 1964 to 2018 as collected by Harvey and Liu (2019). The data used here was downloaded in March 2021.   
## Harvey and Liu's original data can be accessed through https://docs.google.com/forms/d/e/1FAIpQLSflDTiqr2bFsP5tqW7XStIp0-ikwHqK94dbCoMshLfesvKL2g/viewform?vc=0&c=0&w=1
##  dataset "mfalpha" contains CRSP fund id (wficn),  carhart-four-factor-adjusted fund alphas,  alpha t-values and p-values.

require(data.table)
require(mixtools)  # EM algorithm for mixnormal dist.  by G. McLanchlan
require(qvalue)   #  pi0est function
source("FDR_FUN_github.R")


###########################################################################################################
################# factor zoo:  p-hacking ##########
###########################################################################################################


#################### Figure 4 of the paper: P-hacking
pdf("phacking-tstat.pdf", width = 12, height = 6 )
par(mfrow=c(1,3))  #, mar=c(10, 4.5, 0, 1))
hist(abs(factortv$tvalue)[abs(factortv$tvalue) >= 1.96], n=20, xlab = "Cutoff = 1.96", main="", xaxt="n")
axis(side=1, at=c(2,5,10, 15, 20))
hist(abs(factortv$tvalue)[abs(factortv$tvalue) >= 2], n=20,  xlab = "Cutoff = 2", main="", xaxt="n")
axis(side=1, at=c(2,5,10, 15, 20))
hist(abs(factortv$tvalue)[abs(factortv$tvalue) >= 2.57], n=20,  xlab = "Cutoff = 2.57", main="", xaxt="n", ylim=c(0, 120))
axis(side=1, at=c(2,5,10, 15, 20))
#title("P-hacking",  line = -2, outer = TRUE)
dev.off()


### convert tvalue to zscore : truncated at 2 and 7
truncated.tv = factortv$tvalue[abs(factortv$tvalue)>= 2 & abs(factortv$tvalue) < 7]   
pv =  2*(1-pnorm(abs(truncated.tv)))
ZSCORE <- qnorm(1 - pv)

z = ZSCORE

### using moment match to estimate the mean and standard deviation of the distribution under non-null
library(optimx)
a = min(z)
b = max(z)
pp = c(seq(0.4, 0.95, by=0.05))   # prior belief of pi_0
phest <- c()
for(p0 in pp)  phest = rbind(phest, c(p0, twomomentmatch(p0=p0, a=a,b=b,z=z)))

### Table 7:  t-hurdles under FDR control of 0.01 and 0.05
phackingtab(phest, x=z, FDR=0.01)
phackingtab(phest, x=z, FDR=0.05)
### pi'_0:  the zero alpha proportion in the truncated data.  
newp0 = c()
for(i in 1:length(pp)) newp0= c(newp0, tr.loglik.mixnormal(phest[i,1], mu0=0, sigma0=1, mu1=phest[i,2], sigma1=phest[i,3], x=z, a=a,b=b)$newpi0)
cbind(pp, newp0)



###########################################################################################################
#################         Mutual fund performance      ##########
###########################################################################################################

ZSCORE <- qnorm(1 - mfalpha$pv)
if(any(is.infinite(ZSCORE))) { 
  ZSCORE[ZSCORE == Inf] = max(ZSCORE[!is.infinite(ZSCORE)])*1.01
  ZSCORE[ZSCORE == -Inf] = min(ZSCORE[!is.infinite(ZSCORE)])*1.01
}


#################### Figure 5 of the paper: t-stat for mutual fund returns
pdf("mf_tstat_zscore.pdf", width = 10, height = 6 )
par(mfrow=c(1,2))  #, mar=c(10, 4.5, 0, 1))
hist(mfalpha$tv, n=20, xlab = "t-statistic", main="", ylim=c(0, 500))
hist(ZSCORE, n=20, xlab = "z-score", main="")
dev.off()


## model fitting

## Table 8:  Panel A

ini = pi0est(mfalpha$pv, lambda=0.6)$pi0   # set initial value for p_0 as Storey (2002) method with tuning parmaeter = 0.6
#theoretical null model
out0 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,1), sigma=c(1,2),  mean.constr=c(0, NA), sd.constr=c(1, NA), epsilon = 1e-08, maxit=2000, maxrestarts=40)
summary(out0)
#empirical null model
out1 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,1), sigma=c(1,2), epsilon = 1e-08, maxit=2000, maxrestarts=40,  arbvar=T)
summary(out1)

## BIC:  k*ln(n) - 2ln(L): select the model with the minimum BIC value
BIC0 =  3*log(length(ZSCORE))-2*out0$loglik   
BIC1 =  5*log(length(ZSCORE))-2*out1$loglik
diff(c(BIC0, BIC1))
c(BIC0, BIC1)

## AIC:  k*2 - 2ln(L) : select one with minimum AIC value
AIC0 =  3*2-2*out0$loglik   
AIC1 =  5*2-2*out1$loglik
c(AIC0, AIC1)


## Table 8:  Panel B
out = out0
tab <- matrix(NA, nrow = 5, ncol = 6)
tab[,1] <- c(0.01, 0.05, 0.1, 0.15, 0.2)
for(i in 1:5)
{
  lamb <- targetFDR(FDR=tab[i,1], mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1])   # lambda which is associated with a given FDR level
  k = BE(lambda=lamb, w=1/2, mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1], lossonly=F)
  tab[i, 2:6] <- c(k$FNR, k$thurdle, sum(mfalpha$tv >= k$thurdle)/nrow(mfalpha), sum(mfalpha$tv <= -k$thurdle)/nrow(mfalpha),  
                   sum(abs(mfalpha$tv) >= k$thurdle)/nrow(mfalpha))
}
tab


## Table 8:  Panel C
out = out1
tab <- matrix(NA, nrow = 5, ncol = 6)
tab[,1] <- c(0.01, 0.05, 0.1, 0.15, 0.2)
for(i in 1:5)
{
  lamb <- targetFDR(FDR=tab[i,1], mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1])   # lambda which is associated with a given FDR level
  k = BE(lambda=lamb, w=1/2, mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1], lossonly=F)
  tab[i, 2:6] <- c(k$FNR, k$thurdle, sum(mfalpha$tv >= k$thurdle)/nrow(mfalpha), sum(mfalpha$tv <= -k$thurdle)/nrow(mfalpha),  
                   sum(abs(mfalpha$tv) >= k$thurdle)/nrow(mfalpha))
}
tab


#################### Figure 6: mixture distribution plot for mutual fund returns
library("ggplot2")
library("dplyr")
library("ggpubr")  # to place multiple ggplots on one page

pdf("MF_two_componentmix.pdf", width = 10, height = 6 )
out = out0
m1 <- data.frame(x = out$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins =20,  colour = "black", 
                 fill = "white") +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[1], out$sigma[1], pi = out$lambda[1]),
                colour = "red", lwd = 1.1) +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[2], out$sigma[2], pi = out$lambda[2]),
                colour = "blue", lwd = 1.1) + xlab("z-score") +
  ylab("Density")

out = out1
m2 <- data.frame(x = out$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins =20, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[1], out$sigma[1], pi = out$lambda[1]),
                colour = "red", lwd = 1.1) +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[2], out$sigma[2], pi = out$lambda[2]),
                colour = "blue", lwd = 1.1) + xlab("z-score") +
  ylab("Density")
ggarrange(m1,m2, ncol = 2, nrow = 1)
dev.off()




############ ############ ############ ############ ############ ############ ############ ############ ############ 

### The code below is to produce the results related to Yan & Zhang (2017) 18,000 anomalies.  Please contact Yan & Zhang for their data

############    Yan&Zheng 2017 18000 anomalies

## anomaly.vw:  long-short portfolio under value weighting, the associated p_values and t_values Yan&Zheng 2017
## anomaly.ew:  strategy alphas under equal weighting and the associated p_values and t_values 


#################### Figure 2 of the paper: t-stat for anomaly returns
pdf("an_tstat_ewvw.pdf", width = 10, height = 6 )
par(mfrow=c(1,2))  #, mar=c(10, 4.5, 0, 1))
## Equal-weight
hist(anomaly.ew$tstat, n=20, xlab = "Equal weighting", main="", xaxt="n")
axis(side=1, at=c(-8, -5, -2, 0, 2,5,8))
## value-weight
hist(anomaly.vw$tstat, n=20, xlab = "Value weighting", main="", xaxt="n")
axis(side=1, at=c(-5, -2, 0, 2,5))
dev.off()


############################################################# FDR part:  equal weighted
anomaly = anomaly.ew
pv = anomaly$pv
ZSCORE <- qnorm(1 - pv)
if(any(is.infinite(ZSCORE))) { 
  ZSCORE[ZSCORE == Inf] = max(ZSCORE[!is.infinite(ZSCORE)])*1.01
  ZSCORE[ZSCORE == -Inf] = min(ZSCORE[!is.infinite(ZSCORE)])*1.01
}


### Table 5: Panel A
#### Storey (2002) methods to estimate pi0 as initial value to EM algorithm
ini = pi0est(pv, lambda=0.8)$pi0 
ewout0 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,2), sigma=c(1,2),  mean.constr=c(0, NA), sd.constr=c(1, NA), epsilon = 1e-08, maxit=2000, maxrestarts=40)
summary(ewout0)
ewout1 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,2), sigma=c(1,2), epsilon = 1e-08, maxit=2000, maxrestarts=40,  arbvar=T)
summary(ewout1)

## BIC:  k*ln(n) - 2ln(L): select one with minimum BIC value
BIC0 =  3*log(length(ZSCORE))-2*ewout0$loglik   
BIC1 =  5*log(length(ZSCORE))-2*ewout1$loglik
c(BIC0, BIC1)
diff(c(BIC0, BIC1))

## AIC:  k*2 - 2ln(L) : select one with minimum AIC value
AIC0 =  3*2-2*ewout0$loglik   
AIC1 =  5*2-2*ewout1$loglik
c(AIC0, AIC1)


### Table 5: Panel B
out = ewout1
anomaly = anomaly.ew

tab <- matrix(NA, nrow = 5, ncol = 6)
tab[,1] <- c(0.01, 0.05, 0.1, 0.15, 0.2)
for(i in 1:5)
{
  lamb <- targetFDR(FDR=tab[i,1], mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1])   # lambda which is associated with a given FDR level
  k = BE(lambda=lamb, w=1/2, mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1], lossonly=F)
  tab[i, 2:6] <- c(k$FNR, k$thurdle, sum(anomaly$tstat >= k$thurdle)/nrow(anomaly), sum(anomaly$tstat <= -k$thurdle)/nrow(anomaly),  
                   sum(abs(anomaly$tstat) >= k$thurdle)/nrow(anomaly))
}
tab

### Table 5: Panel C
kappaseq <- c(0.05, seq(0.1,0.6,0.1))  ## different rejection rule
tx <- out$posterior[,1]
k = data.frame(FDRTable(tx,kappaseq))
k[, c(1,2,3,4,5,8)]



############################################################# FDR part:  value weighted
anomaly = anomaly.vw
pv = anomaly$pv
ZSCORE <- qnorm(1 - pv)

### Table 6: Panel A
ini = pi0est(pv, lambda=0.8)$pi0  
vwout0 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,1), sigma=c(1,1),  mean.constr=c(0, NA), sd.constr=c(1, NA), epsilon = 1e-06, maxit=2000, maxrestarts=40)
summary(vwout0)
vwout1 = normalmixEM(ZSCORE, lambda=c(ini, 1-ini),  mu=c(0,1), sigma=c(1,1), epsilon = 1e-06, maxit=10000, maxrestarts=40,  arbvar=T) # for vw ret
summary(vwout1)

## BIC:  k*ln(n) - 2ln(L): select one with minimum BIC value
BIC0 =  3*log(length(ZSCORE))-2*vwout0$loglik   
BIC1 =  5*log(length(ZSCORE))-2*vwout1$loglik
c(BIC0, BIC1)
diff(c(BIC0, BIC1))

## AIC:  k*2 - 2ln(L) : select one with minimum AIC value
AIC0 =  3*2-2*vwout0$loglik   
AIC1 =  5*2-2*vwout1$loglik
c(AIC0, AIC1)

### Table 6: Panel B
out = vwout0
anomaly = anomaly.vw

tab <- matrix(NA, nrow = 5, ncol = 6)
tab[,1] <- c(0.01, 0.05, 0.1, 0.15, 0.2)
for(i in 1:5)
{
  lamb <- targetFDR(FDR=tab[i,1], mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1])   # lambda which is associated with a given FDR level
  k = BE(lambda=lamb, w=1/2, mu0=out$mu[1], s0=out$sigma[1],  mu1=out$mu[2], s1=out$sigma[2], pi0=out$lambda[1], lossonly=F)
  tab[i, 2:6] <- c(k$FNR, k$thurdle, sum(anomaly$tstat >= k$thurdle)/nrow(anomaly), sum(anomaly$tstat <= -k$thurdle)/nrow(anomaly),  
                   sum(abs(anomaly$tstat) >= k$thurdle)/nrow(anomaly))
}
tab

### Table 6: Panel C
kappaseq <- c(0.05, seq(0.1,0.6,0.1))  ## different rejection rule
tx <- out$posterior[,1]
k = data.frame(FDRTable(tx,kappaseq))
k[, c(1,2,3,4,5,8)]


##### Figure 3

library("ggplot2")
library("dplyr")
library("ggpubr")  # to place multiple ggplots on one page

pdf("AN_two_componentmix_both.pdf", width = 10, height = 6 )

out = vwout0
m1 <- data.frame(x = out$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins =20,  colour = "black", 
                 fill = "white") +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[1], out$sigma[1], pi = out$lambda[1]),
                colour = "red", lwd = 1.1) +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[2], out$sigma[2], pi = out$lambda[2]),
                colour = "blue", lwd = 1.1) + xlab("z-score") +
  ylab("Density")

out = ewout1
m2 <- data.frame(x = out$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), bins =20, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[1], out$sigma[1], pi = out$lambda[1]),
                colour = "red", lwd = 1.1) +
  stat_function(geom = "line",  fun = plot_mix_comps,
                args = list(out$mu[2], out$sigma[2], pi = out$lambda[2]),
                colour = "blue", lwd = 1.1) + xlab("z-score") +
  ylab("Density")

ggarrange(m2,m1, ncol = 2, nrow = 1)
dev.off()








