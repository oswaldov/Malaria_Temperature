######################################################
## Code to reproduce the consistency analysis 
## requires datafile "data/Prevalence_Pfalciparum.csv"
## code by Oswaldo C. Villena, Leah R. Johnson
## 11 January 2022
######################################################

## Example to validate the model for An. gambiae/P. falciparum 
## First load the necessary R packages and files for the computation and
## visualization 
library("statmod")
library("visreg")
library("dplyr")
library("sjPlot")
library("texreg")
library("coefplot")

## Load data
dat1<-read.csv("data/Prevalence_Pfalciparum.csv", header=TRUE)

## select data >= 1990
datf<-dat1[dat1$year_start>=1990,]

## Select covariates of interest

dat2 <- datf[,c(2,4,5,12,13,14,15,17,18,19)]
summary(dat2)

## Remove NA's 

ww<- which(!is.na(dat2$GDP_adjusted))
dat2<- dat2[ww,]

ww<- which(!is.na(dat2$lGDP))
dat2<- dat2[ww,]


## Select data only for Africa
dat2$AFR <- dat2$continent=="Africa"


## We consider a few combinations of S(T) measures and population
dat2$R0G0.lpop<-dat2$STGZ*log(dat2$pop_den)
dat2$R0G0.lgdp<- dat2$STGZ*log(dat2$GDP_adjusted)
dat2$lR0G0.pop<-log((dat2$STGZ*dat2$pop_den)+1)
dat2$lR0G0.gdp<- log((dat2$STGZ*dat2$GDP_adjusted)+1)

##Check correlation
vary<- c(2,3,7:15)
pairs(dat2[,vary])

## Models for the Presence/Absence Data
## 1 lpop_den only
reg1<-glm(pa ~ lpop_den*AFR - AFR, family=binomial, data=dat2)
summary(reg1)

## 2 lGDP only
reg2<-glm(pa ~ lGDP*AFR - AFR, family=binomial, data=dat2)
summary(reg2)

## 3 ROGO only
reg3<-glm(pa ~ STGZ*AFR - AFR, family=binomial, data=dat2)
summary(reg3)

## 4 socio factors only
reg4<-glm(pa ~ (lGDP +lpop_den)*AFR - AFR, family=binomial, data=dat2)
summary(reg4)

## 5 ROGO + lGDP + lpop_den
reg5<-glm(pa ~ (STGZ+lGDP +lpop_den)*AFR - AFR, family=binomial, data=dat2)
summary(reg5)

## 6 log(ROGO*lpop_den)
reg6<-glm(pa ~ lR0G0.pop*AFR - AFR, family=binomial, data=dat2)
summary(reg6)

## 7 log(ROGO*gdp)
reg7<-glm(pa ~ lR0G0.gdp*AFR - AFR, family=binomial, data=dat2)
summary(reg7)


## BIC - Model comparisons
bics <- c(BIC(reg1), BIC(reg2), BIC(reg3), BIC(reg4), BIC(reg5), BIC(reg6), BIC(reg7)) 


#### Model Comparisons (within sample)

bics
probs <- round(ebics/sum(ebics), 5)
probs


## Proportion of Deviance  explained, compared to null model
D2test<- c(Dsquared(reg1), Dsquared(reg2), Dsquared(reg3), Dsquared(reg4),
           Dsquared(reg5), Dsquared(reg6), Dsquared(reg7))
round(D2test,3)


## Plots residual diagnostics
models<-list(reg1=reg1, reg2=reg2, reg3=reg3, reg4=reg4, reg5=reg5, reg6=reg6, reg7=reg7)

## randomized quantile residuals

qr5<-qresiduals(reg5)
length(qr)
rf<- fitted.values(reg5)
length(rf)
par(mfrow=c(1,2))
plot(rf, qr5, col=as.numeric(dat2$pa)+1, ylim=c(-4,4), xlab="fitted",
     ylab="Randomized quantile residuals", main = expression(atop('Prevalence '*italic(P.~falciparum), ' Africa - Asia')),
     cex.main=1.2, font.main=1)
abline(h=0)
legend(-0.3, 4.7, "a", cex = 2.5, xpd=TRUE, bty="n")
qqnorm(qr, main="Normal Q-Q Plot \nAfrica - Asia", ylim=c(-4,4), cex.main=1.2, font.main=1)
abline(0,1, col=2)
legend(-6, 4.7, "b", cex=2.5, xpd=TRUE, bty="n")



### Visualizations of the model fit.
## The posterior probability that S(T)>0 vs the probability of infection 
## caused by P. falciparum and transmited by An. gambiae


par(mfrow=c(1,1))
visreg(reg5, "STGZ", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (S(T)>0)), ylab=expression(Probability~of~infection~italic(P.~falciparum)), cex.lab=1.1,
       main=expression(Prevalence~italic(Plasmodium~falciparum)), legend=FALSE)
legend(0.7, 0.15, 
       c("AFRICA", "ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend(-0.18, 1.05, paste('',"a"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)

## The natural log of the probability of S(T)>0 * the percapita GDP vs the probability of infection 
## caused by P. falciparum and transmited by An. gambiae

visreg(reg7, "lR0G0.gdp", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(log(S(T)>0%*%GDP)), ylab=expression(probability~of~infection~italic(P.~falciparum)), cex.lab=1.1,
       main=expression(Prevalence~italic(Plasmodium~falciparum)), legend=FALSE)
legend(7, 0.15, 
       c("AFRICA", "ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend(-1.7, 1.04, paste('',"b"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)

## Multiple plots of the probability s(T)>0 vs the probability of infection 
## caused by P. falciparum and transmited by An. gambiae

my.GDP<-c(min(dat2$lGDP), median(dat2$lGDP), max(dat2$lGDP))
my.GDP

my.lpop<-c(min(dat2$lpop_den), median((dat2$lpop_den)-1), max((dat2$lpop_den)-1))
my.lpop


## Plot

par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "STGZ", scale="response", by="AFR", rug=2, overlay=TRUE, 
           cond=list(lpop_den=my.lpop[i], lGDP=my.GDP[j]), 
           ylim=c(0,1),line=list(col=c("red","dodgerblue"), lwd=3),
           fill=list(col=c("#FF4E3780","#37b6ff")), legend=FALSE,
           xlab=expression(Prob (S(T)>0)), ylab=expression(Probability~of~infection~italic(P.~falciparum)), cex.lab=0.9)
    legend("topleft",
           c(paste("pop density = ", round(10^(my.lpop[i])), "K", sep=""),
             paste("GDP =", round(10^my.GDP[j]), " USD", sep ="")),
           bty="n", y.intersp=1)
    legend("bottomright", 
           c("AFR","ASIA"), cex = 1.2, text.font =3,
           col=c("dodgerblue","red"), lwd=c(3,3),
           bty="n", y.intersp=1)
  }
}


## Posterior samples differences between mosquito/parasite systems
## Here we showed the code to calculate the sample differences of the S(T) posterior distributions at the minimum thermal limit

## Load the data that contains the minimum, optimum, and maximum values from the S(T) metric


aspf<- load("data/STPosterior_AnSteph_Pfalc.Rsave")

aspv<- load("data/STPosterior_AnSteph_Pvivax.Rsave")

agpf<- load("data/STPosterior_AnGamb_Pfalc.Rsave")

agpv<- load("data/STPosterior_AnGamb_Pvivax.Rsave")

## Comparison for An. stephensi/P. falciaprum vs An. stephensi/P. vivax 

vec1<- R0minASPF; vec2<- R0minASPV

N<- length(vec1) ## length of the vector

cat <- function(x,y){  ## function that calculates the difference between compared vectors
  diffe <- x - y
  return(diffe)
}

alst<- cat(vec1,vec2) ## Apply cat function to obtain values from the vectors difference. 


## percentages
countgre<- (sum(alst > 0)/N) ## counts greater than zero

countless<- (sum(alst < 0)/N) ## counts less than zero

counteq<- (sum(alst == 0)/N) ## counts equal to zero


## Comparison for An. gambiae/P. falciaprum vs An. gambiae/P. vivax 

vec1<- R0minAGPF; vec2<- R0minAGPV

algb<- cat(vec1,vec2)


## percentages
countgre<- (sum(algb > 0)/N)

countless<- (sum(algb < 0)/N)

counteq<- (sum( algb == 0)/N)


## Comparison for An. stephensi/P. falciparum vs An. gambiae/P. falciparum 

vec1<- R0minASPF; vec2<- R0minAGPV

alsfgf<- cat(vec1,vec2)

## percentages
countgre<- (sum(alsfgf > 0)/N)

countless<- (sum(alsfgf < 0)/N)

counteq<- (sum( alsfgf == 0)/N)


## Comparison for An. stephensi/P. vivax vs An. gambiae/P. vivax 

vec1<- R0minASPV; vec2<- R0minAGPV

alsvgv<- cat(vec1,vec2)

## percentages
countgre<- (sum(alsvgv > 0)/N)

countless<- (sum(alsvgv < 0)/N)

counteq<- (sum( alsvgv == 0)/N)

## Plot the four graphs for the comparison of S(T) posterior distributions at the minimum thermal limit

op<- par(
  oma=c(0,0,2,0),
  mfrow = c(2,2))
plot(alst, ylim=c(-10,8),
     xlab="S(T) posterior distributions",ylab="Temperature differences (°C)", cex.lab=1.2,
     col = ifelse(alst > 0,'blue','red'))
legend(-1000, 9.5,("a"),bty='n', cex = 1.7)
abline(h=0, col="black", lwd=2)
legend("topright", inset=c(0.01,0.83), legend=c("positive","negative"), col=c("blue","red"),
       pch=c(1,1), title="An. steph/Pfalc - An. Steph/Pviv", cex=0.8, horiz = TRUE)


plot(algb, ylim=c(-10,8), 
     xlab="S(T) posterior distributions",ylab="Temperature differences (°C)", cex.lab=1.2,
     col = ifelse(algb > 0,'blue','red'))
legend(-1000, 9.5,("b"),bty='n', cex = 1.7)
abline(h=0, col="black", lwd=2)
legend("topright", inset=c(0.01,0.83), legend=c("positive","negative"), col=c("blue","red"),
       pch=c(1,1), title="An. gamb/Pfalc - An. gamb/Pviv", cex=0.8, horiz = TRUE)

plot(alsfgf, ylim=c(-10,8), 
     xlab="S(T) posterior distributions",ylab="Temperature differences (°C)", cex.lab=1.2,
     col = ifelse(alsfgf > 0,'blue','red'))
legend(-1000, 9.5,("c"),bty='n', cex = 1.7)
abline(h=0, col="black", lwd=2)
legend("topright", inset=c(0.01,0.83), legend=c("positive","negative"), col=c("blue","red"),
       pch=c(1,1), title="An. steph/Pfalc - An. gamb/Pfalc", cex=0.8, horiz = TRUE)

plot(alsvgv, ylim=c(-10,8), 
     xlab="S(T) posterior distribution",ylab="Temperature differences (°C)", cex.lab=1.2,
     col = ifelse(alsvgv > 0,'blue','red'))
legend(-1000, 9.5,("d"),bty='n', cex = 1.7)
abline(h=0, col="black", lwd=2)
legend("topright", inset=c(0.01,0.83), legend=c("positive","negative"), col=c("blue","red"),
       pch=c(1,1), title="An. steph/Pviv - An. gamb/Pviv", cex=0.8, horiz = TRUE)

par(op)
op<- par(usr=c(0,1,0,1),
         xpd=NA)
mtext("Comparison of S(T) posterior distributions - Minimum temperature limit     ", side=3, line= 1, cex=1)
