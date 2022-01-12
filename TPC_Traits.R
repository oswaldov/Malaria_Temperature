###################################################################
## Code to reproduce the estimation of the thermal performance curves (TPC)
## requires datafile "data/traits.csv"
## code by Oswaldo C. Villena, Leah R. Johnson
## 11 January 2022
###################################################################


## First load the necessary R packages and files for the computation and
## visualization
library('rjags')
library('MASS')


## Next load data
data.all<-read.csv("data/traits.csv", header=TRUE)

## Select the trait of interest. Here we showed as an example MDR for An. gambiae
data<-data.all[which(data.all$trait.name=="mdr" & data.all$specie=="An. gambiae"), 2:9]
plot(data$T, data$trait) ## visualize your data


## specify the parameters that control the MCMC
n.chains<-5
n.adapt<-10000
n.samps<-20000
n.burn<- 10000


##Choose the appropiate model for the specific trait 

##Briere model (MDR, PDR, a)

jags_briere.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)T(0,)
mu.temp[i] <- c*T[i]*(T[i]-T0)*sqrt((Tm-T[i])*(Tm>T[i]))
mu[i] <- 0*(mu.temp[i]<0) + mu.temp[i]*(mu.temp[i]>0)

}

c ~ dgamma(1,10)
Tm ~ dunif(25,45)
T0  ~ dunif(0, 24)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"


## Concave down quadratic model (PEA, EFD, bc)

jags_quad.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)T(0,)
mu[i] <- -qd*(T[i]-T0)*(T[i]-Tm)*((T[i]>T0))*((T[i]<Tm))
}

Tm  ~ dunif(25,45)
T0 ~ dunif(0,24)
qd  ~ dgamma(1,1)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"


## Concave up quadratic model (mu)

jags_quad.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)T(0,)
mu[i] <- inter-n.slope*T[i]+qd*T[i]^2
}

inter ~ dgamma(2,2)
n.slope ~ dgamma(1, 1)
qd  ~ dgamma(2,2)
sigma<-1/tau
tau ~ dnorm(1000, 1/500)

}"



## Use jags.model for the specific model with the appropiate
## default priors
jags1 <- jags.model(textConnection(jags_briere.bug),  ## change the model according to the trait. This example is for MDR
                    data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                    n.chains = n.chains, inits=list(Tm=31, T0=5, c=0.00007),
                    n.adapt = n.adapt) 
update(jags1, n.burn)


## The coda.samples() functions takes n.samps new samples, and saves
## them in the coda format, which we use for visualization and
## analysis

coda.samps1 <- coda.samples(jags1, c('c','Tm', 'T0', 'sigma'), n.samps)
summary(coda.samps1)

## check for convergence
plot(coda.samps1,ask=T)


## This command combines the samples from the n.chains into a format
## that we can use for further analyses. Use appropiate model for specific traits
samps<-make.briere.samps(coda.samps1, nchains=n.chains, samp.lims=c(1, n.samps))
head(samps)
samps$tau<-1/samps$sigma


## Next we want to use the parameter samples to get posterior samples
## of the temperature rsponses themselves
Temps<-seq(0,50, by=0.1)
out<-make.sims.temp.resp(sim="briere", MDR.sampsgampf, Temps, thinned=seq(1,n.samps, length = 1000)) ## Example for MDR
summary(out)

## and then we calculate the 95% inner quantile/HPD
q<-temp.sim.quants(out$fits, length(Temps))



## We can then plot the data with the fits/quantiles

par(mfrow=c(1,1),
    mar = c(5.1,4.8,4.1,2.1))
plot(data$T, data$trait, xlim = c(0, 45), ylim = c(0, 0.15), cex.lab=1.6, cex.axis=1.3,
     pch=19, cex=1.5,
     xlab="Temperature (C)",
     ylab="Mosquito development rate, MDR",
     box(lty = "solid", bty = "o"))

legend(1,0.16, paste('', "f"), text.font = 3, bty = "n", xjust= 0, cex= 2.5)

add.sim.lines(Temps, sim.data=out$fits, mycol=1, lwd=2)
lines(Temps, hpdl, col="blue", lty=5, lwd=2)
lines(Temps, hpdh, col="blue", lty=5, lwd=2)




