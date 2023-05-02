## name: Ziwen Zhong | student number: s2326022

## Overview
## Code to calculate the excess deaths in the UK from 2020 to 2022.

## We predict the deaths from 2020 to 2022 using death rate for each age class
## and sex based on the previous years, and modify them by multiplying weekly 
## coefficient d, which also come from previous data, out of the consideration 
## for seasonal variation in mortality rates.
## The excess deaths is the difference between the observed deaths and the 
## predicted deaths.
## Then we build a simple time series model to fit the excess deaths, and use
## Gibbs sampling to obtain the samples of 2 parameters of the model that we 
## are interested in.


###############################################################################
# load required packages
library('ggplot2')
library('rjags')


## load data
pop <- read.table('lt1720uk.dat')
death <- read.table('death1722uk.dat')
## attach two tables
attach(pop); attach(death)


death_predictor <- function(Nf,Nm,mf,mm,d) {
  ## summary
  ##   a function to predict total number of deaths weekly
  ## input
  ##   Nf, Nm: integer. population of female and male for each age classes
  ##   mf, mm: mortality rate of female and male for each age classes
  ##   d: mortality rate modifier for each week
  ## output
  ##   a vector of integer. predict deaths for each week
  qf <- 1-exp(-mf/52)   ## calculate deaths proportion for female
  qm <- 1-exp(-mm/52)   ## calculate deaths proportion for male
  pop_death <- rep(0,length(d))   ## a vector to store predicted deaths weekly
  for (j in 1:length(d)) {
    Df <- 0.9885*d[j]*qf*Nf   ## calculate deaths for female
    Dm <- 0.9885*d[j]*qm*Nm   ## calculate deaths for male
    Nf_star <- Nf - Df   ## population of female before ageing process
    Nm_star <- Nm - Dm   ## population of male before ageing process
    ## apply the ageing process, 1/52 people from the former age class move in
    ## and 1/52 people in this age class move out
    Nf[-1] <- 51/52*Nf_star[-1] + 1/52*Nf_star[-length(Nf_star)]
    Nm[-1] <- 51/52*Nm_star[-1] + 1/52*Nm_star[-length(Nm_star)]
    pop_death[j] <- sum(Df+Dm)  ## record deaths
  }
  pop_death
}


## calculate the overall excess deaths start from 2020 and within 2020
pred_deaths <- death_predictor(fpop20,mpop20,mf,mm,d[157:length(d)])
excess_deaths <- sum(deaths[157:length(d)]) - sum(pred_deaths)
excess_deaths_2020 <- sum(deaths[157:208]) - sum(pred_deaths[1:52])


## 1st plot: Observed and Predicted Deaths against Week
week_from_2020 <- 157:length(d)-156     ## reset week index 
plot(week_from_2020, deaths[157:length(d)],   ## plot observed deaths
     xlab='Week from 2020',ylab='Deaths',
     ## report the excess deaths in title
     main=paste('Excess deaths in 2020 is',round(excess_deaths_2020,2),
                ', overall is',round(excess_deaths,2)))    
lines(week_from_2020,pred_deaths, col=2)     ## add curve of predict deaths
legend(x='topright', legend=c('observed deaths','predicted deaths'),
       fill=c(1,2),inset=.02,cex=0.8)     ## add the legend


## 2nd plot: Cumulative Excess Deaths against Week
## calculate excess deaths for each week
excess_deaths_weekly <- deaths[157:length(d)] - pred_deaths
## calculate cumulative excess deaths fot each week
cumn_excess_deaths <- cumsum(excess_deaths_weekly)
plot(week_from_2020,cumn_excess_deaths,
     xlab='Week from 2020',ylab='Cumulative excess deaths',
     main='Cumulative excess deaths against Week')


## use a MCMC model defined in 'model.jags' to fit the weekly excess deaths
x = excess_deaths_weekly
x[c(51,52,53,105,106)] <- NA  ## set abnormal data to NA (unused to fit model)
mod <- jags.model('model.jags',data=list(x=x,N=length(x)))  ## fit model
## sample 10000 mu and rho from the model
sam <- jags.samples(mod,c("mu","rho"),n.iter=10000)


## 3rd & 4th plot: the trace plot and the histograms of rho
rho <- sam$rho     ## extract rho
## extract the start and the end of iterations of rho
iterations <- attr(sam$rho,"iterations")
iterations <- seq(iterations['start'],iterations['end'],iterations['thin'])
## plot the trace plot of rho
plot(iterations,rho,type='l',xlab='Iterations',ylab='',
     main=expression(paste('Trace of ',rho)))
## plot the histograms of rho
hist(rho,main=expression(paste('Histogram of ',rho)),xlab=expression(rho))


## extract mu and calculate the expected mu
mu <- matrix(sam$mu,nrow(sam$mu),ncol(sam$mu))
mu_expect <- apply(mu,1,mean)


## 5th plot: Sampled mu, Expected mu, and Excess Deaths against Week
## empty plot
plot(week_from_2020,mu_expect,type='n',xlab='Week from 2020',
     ylab='Excess deaths',main='Excess deaths against Week')
## add curves of sampled mu (every 50th)
for (j in seq(50, 10000, 50)) lines(week_from_2020,mu[,j],col=8)
## add curve of expected of mu
lines(week_from_2020,mu_expect,type='l',col=4)
## add excess death point (used)
points(week_from_2020[-c(51,52,53,105,106)],
       excess_deaths_weekly[-c(51,52,53,105,106)])
## add excess death point (unused) in red colour
points(week_from_2020[c(51,52,53,105,106)],
       excess_deaths_weekly[c(51,52,53,105,106)],col=2)
## add legend
legend(x='topright',legend=c(expression(paste('sampled ',mu)),
                             expression(E(mu)),
                             expression(paste('observed excess deaths ',x[i])),
                             expression(paste('unused ',x[i]))),
       fill=c(8,4,1,2),inset=.02,cex=0.8)


## 6th plot: Residual (x - Expected mu) against Week
plot(week_from_2020,x-mu_expect,type='p',xlab='Week from 2020',
     ylab='Residuals',main='Residuals against Week')

## detach two tables
detach(pop); detach(death)









