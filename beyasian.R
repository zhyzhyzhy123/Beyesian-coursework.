library(ggplot2)
library(ggpubr)
library(patchwork)
setwd("~/Desktop/modules/Bayesian/cw")
#need to change path if run on your own computer
X <- read.csv('cw.csv')[,2]
##############Data Generation
N <- 100000
##############Posterior 1
res1 <- matrix(0,N,2)
res1[1,] <- c(3600, 250)
for (i in 2:N) {
  #I choose folded normal distribution as proposals for both parameters
  mustar <- abs(rnorm(1, res1[i-1,1],25))
  sigmastar <- abs(rnorm(1, res1[i-1,2],25))
  #loglikelihood
  ratiotop <- sum(dnorm(X, mustar, sigmastar, log = T))
  #prior
  ratiotop <- ratiotop + dbeta((mustar-3000)/600, 4.86, 3.44, log = T)+ dbeta((sigmastar-50)/200, 3, 3, log = T)
  #loglikeihood
  ratiobot <- sum(dnorm(X, res1[i-1,1], res1[i-1,2], log = T))
  #prior
  ratiobot <- ratiobot + dbeta((res1[i-1,1]-3000)/600, 4.86, 3.44, log = T)+ dbeta((res1[i-1,2]-50)/200, 3, 3, log = T)
  if(mustar>=3000 & mustar<=3600 & sigmastar <= 250 & sigmastar >= 50){
    u <- runif(1)
    if(log(u)<(ratiotop-ratiobot)){
      res1[i,] <- c(mustar, sigmastar)
    }else{
      res1[i,] <- res1[i-1,]
    }
  }else{
    res1[i,] <- res1[i-1,]
  }
}
###############Posterior 2
res2 <- matrix(0,N,2)
res2[1,] <- c(3600, 250)
for (i in 2:N) {
  #I choose folded normal distribution as proposals for both parameters as well
  mustar <- abs(rnorm(1, res2[i-1,1],25))
  sigmastar <- abs(rnorm(1, res2[i-1,2],25))
  #loglikelihood
  ratiotop <- sum(dnorm(X, mustar, sigmastar, log = T))
  #prior
  #ratiotop <- ratiotop + dbeta((mustar-3000)/600, 4.86, 3.44, log = T)+ dbeta((sigmastar-50)/200, 3, 3, log = T)
  #loglikeihood
  ratiobot <- sum(dnorm(X, res2[i-1,1], res2[i-1,2], log = T))
  #prior
  #ratiobot <- ratiobot + dbeta((res1[i-1,1]-3000)/600, 4.86, 3.44, log = T)+ dbeta((res1[i-1,2]-50)/200, 3, 3, log = T)
  if(mustar>=3000 & mustar<=3600 & sigmastar <= 250 & sigmastar >= 50){
    u <- runif(1)
    if(log(u)<(ratiotop-ratiobot)){
      res2[i,] <- c(mustar, sigmastar)
    }else{
      res2[i,] <- res2[i-1,]
    }
  }else{
    res2[i,] <- res2[i-1,]
  }
}
###########Analysics for posteriot 1
#drop the burn-in period in 0-10000
mean(res1[10000:100000,1])
var(res1[10000:100000,1])
mean(res1[10000:100000,2])
var(res1[10000:100000,2])
###########Analysics for posteriot 2
mean(res2[10000:100000,1])
var(res2[10000:100000,1])
mean(res2[10000:100000,2])
var(res2[10000:100000,2])
###########Ploting
setwd("~/Desktop/modules/Bayesian/cw/tex")
df1 <- data.frame(res1,index = 1:100000)
colnames(df1)[1:2] <- c('mu', 'sigma')
df2 <- data.frame(res2,index = 1:100000)
colnames(df2)[1:2] <- c('mu', 'sigma')
##########Posterior 1
p11 <- ggplot(df1, aes(x = index, y = mu))+
  geom_line()+labs(x = 'i', y = expression(mu[1]), title = expression('Line Chart of '*mu[1]))+
  theme(plot.title = element_text(hjust = 0.5))
#line chart of mu1
p12 <- ggAcf(df1[,1])+labs(title = expression('ACF of '*mu[1]))+
  theme(plot.title = element_text(hjust = 0.5))
#acf plot of mu1
p14 <- ggplot(df1, aes(x = index, y = sigma))+
  geom_line()+labs(x = 'i', y = expression(sigma[1]), title = expression('Line Chart of '*sigma[1]))+
  theme(plot.title = element_text(hjust = 0.5))
#line chart of sigma1
p15 <- ggAcf(df1[,2])+labs(title = expression('ACF of '*sigma[1]))+
  theme(plot.title = element_text(hjust = 0.5))
#acf plot of sigma1
p1 <- ggarrange(p11,p12, p14, p15, ncol=2,nrow=2)+
  plot_annotation(title = 'Samples from the First Posterior')&
  theme(plot.title = element_text(hjust = 0.5))
#merge all plot into one
ggsave('f1.png', width = 6.9, height = 4)
###########Posterior 2
p21 <- ggplot(df2, aes(x = index, y = mu))+
  geom_line()+labs(x = 'i', y = expression(mu[2]), title = expression('Line Chart of '*mu[2]))+
  theme(plot.title = element_text(hjust = 0.5))
#line chart of mu2
p22 <- ggAcf(df2[,1])+labs(title = expression('ACF of '*mu[2]))+
  theme(plot.title = element_text(hjust = 0.5))
#acf plot of mu2
p24 <- ggplot(df2, aes(x = index, y = sigma))+
  geom_line()+labs(x = 'i', y = expression(sigma[2]), title = expression('Line Chart of '*sigma[2]))+
  theme(plot.title = element_text(hjust = 0.5))
#line chart of sigma2
p25 <- ggAcf(df2[,2])+labs(title = expression('ACF of '*sigma[2]))+
  theme(plot.title = element_text(hjust = 0.5))
#acf plot of sigma2
p2 <- ggarrange(p21,p22, p24, p25, ncol=2,nrow=2)+
  plot_annotation(title = 'Samples from the First Posterior')&
  theme(plot.title = element_text(hjust = 0.5))
#merge all plots
ggsave('f2.png', width = 6.9, height = 4)
##########################predicted distribution
###########distribution 1
pre1 <- NULL
for (i in 10000:100000) {
  mup <- res1[,1][i]
  sigmap <- res1[,2][i]
  pre1 <- c(pre1, rnorm(1, mup, sigmap))
}
###########distribution 2
pre2 <- NULL
for (i in 10000:100000) {
  mup <- res2[,1][i]
  sigmap <- res2[,2][i]
  pre2 <- c(pre2, rnorm(1, mup, sigmap))
}
############plot
#draw predicted distribution
dfpre <- data.frame(p1 = pre1, p2 = pre2)
p31 <- ggplot(dfpre,aes(x = p1, y = ..density..))+
  geom_histogram(bins = 50, fill = "cornsilk", colour = "grey60")+
  xlab('x')+ylab(expression(f[1](x)))
p32 <- ggplot(dfpre,aes(x = p2, y = ..density..))+
  geom_histogram(bins = 50, fill = "cornsilk", colour = "grey60")+
  xlab('x')+ylab(expression(f[2](x)))
p3 <- ggarrange(p31, p32, nrow=2)+
  plot_annotation(title = 'Predictive Distribution from Two Posteriors')&
  theme(plot.title = element_text(hjust = 0.5))
ggsave('f3.png', width = 6.9, height = 4)
