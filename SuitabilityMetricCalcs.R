###################################################################
## Code to reproduce the estimation of the suitability metric S(T)
## requires datafile "Angamb_Pfalc_samps.Rsave"
## code by Oswaldo C. Villena, Leah R. Johnson
## 11 January 2022
###################################################################


##Example to calculate the  S(T) for An. gambiae/P. falciparum
## First load the necessary R packages and files for the computation and
## visualization 
library(matrixStats)
library(Hmisc)


## First, we identify the files where the posterior samples for each
## component of R0 are saved.

post<- load("data/Angamb_Pfalc_samps.Rsave") ## Rdata file from previous step


## Next we set up the temperatures over which we will be evaluating S(T)
## as well as the thinning interval for the samples.

temp = seq(5,45,by=0.1)	##temperature sequence

t<-length(temp)  ## length of temperature

n = dim(a.sampsgampf)[1]    	## length of samples

thinned<-seq(1, n, by=20)  ## thinned samples

lthin<-length(thinned)  ## number of thinned samples

ec<-0.000001  ## small constant used to keep denominators from being numerically zero


## Function encoding the value of S(T) as a function of the parameters
myST <-function(a, PDR, MDR, efd, pea, bc, mu){
  ((a^2*bc*(efd*pea*MDR/mu^2)*exp(-mu/(PDR+ec)))/(mu+ec))^0.5
}

## The following code runs through the samples and calculates the
## posterior trajectories (i.e. as a function of temperatures) of each
## component and of S(T), and sets up a matrix to store them (each
## column is a trajectory)

ST <-matrix(NA,t,lthin)
a<-PDR<-MDR<-efd<-e2a<-bc<-mu<-matrix(NA,t,lthin)
for (j in 1:lthin){
  if(j%%50==0) cat("iteration", j, "\n")
  ## calculate parameter trajectories
  i<-thinned[j]
  a[,j] = briere(temp,a.sampsgampf[i,3],a.sampsgampf[i,2],a.sampsgampf[i,1])
  PDR[,j] = briere(temp,PDR.sampsgampf[i,3],PDR.sampsgampf[i,2],PDR.sampsgampf[i,1])
  MDR[,j] = briere(temp,MDR.sampsgampf[i,3],MDR.sampsgampf[i,2],MDR.sampsgampf[i,1])
  efd[,j] = quad.2(temp,efd.sampsgampf[i,1],efd.sampsgampf[i,2],efd.sampsgampf[i,3])
  pea[,j] = quad.2.trunc(temp,pea.sampsgampf[i,1],pea.sampsgampf[i,2],pea.sampsgampf[i,3]) 
  bc[,j] = quad.2.trunc(temp,bc.sampsgampf[i,1],bc.sampsgampf[i,2],bc.sampsgampf[i,3])
  mu[,j] = quad(temp,mu.sampsgampf[i,1], -mu.sampsgampf[i,2],mu.sampsgampf[i,3])
  
  
## Calculate S(T)
  ST[,j]<-myST(a[,j], PDR[,j], MDR[,j], efd[,j], pea[,j], bc[,j], mu[,j])
  
}


## Next calculate the posterior mean trajectory of each component of
## S(T). These will be used as part of the uncertainty
## analysis.
a.M<-rowMeans(a)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
efd.M<-rowMeans(efd)
pea.M<-rowMeans(e2a)
bc.M<-rowMeans(bc)
mu.M<-rowMeans(mu)
ST.M<-rowMeans(ST)


## Build matrices to hold results
ST.a<-ST.bc<-ST.efd<-ST.pea<-ST.MDR<-ST.mu<-ST.PDR<-matrix(NA,t,lthin)


## For uncertainty analysis: calculate posterior samples for S(T) with
## all but a single component fixed the posterior mean
for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  ## calculate derivative trajectories
  i<-thinned[j]
  ## Calculate R0 with most components set to their means
  ST.a[,j] = myST(a[,j], PDR.M, MDR.M, efd.M, pea.M, bc.M, mu.M)
  ST.bc[,j] = myST(a.M, PDR.M, MDR.M, efd.M, pea.M, bc[,j], mu.M)
  ST.efd[,j] = myST(a.M, PDR.M, MDR.M, efd[,j], pea.M, bc.M, mu.M)
  ST.pea[,j] = myST(a.M, PDR.M, MDR.M, efd.M, pea[,j], bc.M, mu.M)
  ST.MDR[,j] = myST(a.M, PDR.M, MDR[,j], efd.M, pea.M, bc.M, mu.M)
  ST.mu[,j] =myST(a.M, PDR.M, MDR.M, efd.M, pea.M, bc.M, mu[,j])
  ST.PDR[,j] = myST(a.M, PDR[,j], MDR.M, efd.M, pea.M, bc.M, mu.M)
}


## Calculate the distance within the inner 95% quantile for S(T) overall
## (ST.q) and for the posterior of S(T) with each component held fixed
ST.q<-  apply(ST, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST, 1, FUN=quantile, probs=0.025,na.rm = F)

a.q<-  apply(ST.a, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.a, 1, FUN=quantile, probs=0.025,na.rm = F)
bc.q<- apply(ST.bc, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.bc, 1, FUN=quantile, probs=0.025, na.rm = F)
efd.q<- apply(ST.efd, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.efd, 1, FUN=quantile, probs=0.025, na.rm= F)
pea.q<-apply(ST.pea, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.pea, 1, FUN=quantile, probs=0.025, na.rm = F)
MDR.q<-  apply(ST.MDR, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.MDR, 1, FUN=quantile, probs=0.025, na.rm= F)
mu.q <-  apply(ST.mu, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.mu, 1, FUN=quantile, probs=0.025, na.rm = F)
PDR.q<- apply(ST.PDR, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(ST.PDR, 1, FUN=quantile, probs=0.025, na.rm=F)

## Next plot relative width of quantiles


plot(temp, a.q/(ST.q +ec), col=2, type="l", ylim=c(0,1), lwd=3,
     xlab="Temperature (°C)", ylab="Relative width of HPD intervals", xlim=c(15,38))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
lines(temp, bc.q/(ST.q +ec), col=3, lwd=3)
lines(temp, efd.q/(ST.q +ec), col=4, lwd=3)
lines(temp, pea.q/(ST.q +ec), col=5, lwd=3)
lines(temp, MDR.q/(ST.q +ec), col=6, lwd=3)
lines(temp, mu.q/(ST.q +ec), col=7, lwd=3)
lines(temp, PDR.q/(ST.q +ec), col=8, lwd=3)
## Add legend
leg.text<-c("a", "bc", "EFD", "pea", "MDR", "mu", "PDR")
leg.col<-seq(2, 8, by=1)
legend(15, 1,  leg.text, col=leg.col, lwd=c(3,3,3,3,3,3,3))

legend("topright", paste('', "A"), text.font = 3, bty = "n", cex= 1.8)

## Calculate the distribution of the lower and upper limits of S(T) and peak S(T).
ST.min<-ST.max<-ST.peak<-rep(NA, length(thinned))
## PEAK
for(i in 1:length(thinned)){
  ww<-which(ST[,i]==max(ST[,i]))
  ST.peak[i]<-temp[ww[1]]
}
## MINIMUM
for(i in 1:length(thinned)){
  ww<-which(ST[,i]>0)
  ST.min[i]<-temp[ww[1]-1]
}
## MAXIMUM
for(i in 1:length(thinned)){
  ww<-which(ST[,i]>0)
  lw<-length(ww)
  ST.max[i]<-temp[ww[lw]+1]
}

## plot mean R0 with it's quantiles, all scaled by max mean S(T)
ey <- expression(italic(An.~gambiae/P.~falciparum)~relative~S(T))
par(mar=c(5,5,2,2), mfrow=c(1,1), bty="n")
par(mfrow=c(1,1), bty="n")
ST.scale<-max(ST.M)
ST.q<-temp.sim.quants(ST, length(temp))##/R0.scale
plot(temp, ST.M/ST.scale, type="l", col=1, lwd=3, xlim=c(15, 38), ylim=c(0, 3),
     ylab= ey, xlab="Temperature (°C)", cex.lab=1.2)
box(lty = "solid", bty = "o")
minor.tick(nx=5,ny=5,tick.ratio = 0.5)

add.sim.lines(temp, sim.data=NULL, q=ST.q/ST.scale, mycol=2, lwd = 2)

legend("topright", paste('', "A"), text.font = 3, bty = "n", cex= 1.8)



