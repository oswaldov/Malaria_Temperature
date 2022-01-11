######################################################
## Code to reproduce the consistency analysis 
## requires datafile "data/Prevalence_Pfalciparum.csv"
##code by Oswaldo C. Villena, Leah R. Johnson
##11 January 2022
######################################################

##Example to validate the model for An. gambiae/P. falciparum
## First load the necessary R packages and files for the computation and
## visualization 
library("statmod")
library("visreg")
library("plyr")
library("modEvA")
library(stargazer)
library(sjPlot)
library(texreg)
library(jtools)
library(huxtable)
library(ggstance)
library(broom.mixed)
library(gridExtra)
library(coefplot)
library(dotwhisker)

##load data
dat1<-read.csv("data/Prevalence_Pfalciparum.csv", header=TRUE)

## select data for the appropiate continent and above year 1990
dataAF <- dat1[dat1$continent == "Africa",]

datf<-dataAF[dataAF$year_start>=1990,]
attach(datf)

##Select covariates of interest

myvars<- c("continent", "pop_den", "pf_pos", "pf_pr", "examined", "pf_pa", "gdpvl", "qtmean1", "Roq1", "R0G0qt1")
dat2<- datf[myvars]
summary(dat2)

##remove NA 

ww<- which(!is.na(dat2$gdpvl))
dat2<- dat2[ww,]

ww<- which(!is.na(dat2$pop_den))
dat2<- dat2[ww,]

## log pop density and gdppp

dat2$lpop_den <- log(dat2$pop_den)

dat2$lgdpvl <- log(dat2$gdpvl)

## we consider a few combinations of S(T) with covariates

dat2$l.prev<-log(dat2$pf_pa*dat2$pf_pr+1)
dat2$STGZ.lpop<-dat2$STGZqt1*log(dat2$pop_den)
dat2$STGZ.lgdp<- dat2$STGZqt1*log(dat2$gdpvl)
dat2$lST.pop<- log((dat2$STq1*dat2$pop_den)+1)
dat2$lSTGZ.pop<-log((dat2$STGZqt1*dat2$pop_den)+1)
dat2$lST.gdp<- log((dat2$STq1*dat2$gdpvl)+1)
dat2$lSTGZ.gdp<- log((dat2$STGZqt1*dat2$gdpvl)+1)

## check for correlated variables
vary<- c(2,6,9,10,11,12,17,19)
pairs(dat2[,vary])
round(cor(dat2[,c(2,6,9,10,11,12,17,19)]), digits=2)

## Models for the Presence/Absence Data


## 1 lpop_den only
reg1<-glm(pf_pa ~ lpop_den, family=binomial, data=dat2)
summary(reg1)

## 2 lGDP only
reg2<-glm(pf_pa ~ lgdpvl, family=binomial, data=dat2)
summary(reg2)

## 3 STGZ only
reg3<-glm(pf_pa ~ STGZqt1, family=binomial, data=dat2)
summary(reg3)

## 4 socio-economic factors only
reg4<-glm(pf_pa ~ (lgdpvl +lpop_den), family=binomial, data=dat2)
summary(reg4)

## 5 STGZ + lGDP + lpop_den
reg5<-glm(pf_pa ~ (STGZqt1+lgdpvl+lpop_den), family=binomial, data=dat2)
summary(reg5)

## 6 log(STGZ*lpop_den)
reg6<-glm(pf_pa ~ lSTGZ.pop, family=binomial, data=dat2)
summary(reg6)

## 7 log(STGZ*gdp)
reg7<-glm(pf_pa ~ lSTGZ.gdp, family=binomial, data=dat2)
summary(reg7)

## BIC - Model comparisons
bics <- c(BIC(reg1), BIC(reg2), BIC(reg3), BIC(reg4), BIC(reg5), BIC(reg6), BIC(reg7)) 


#### Model Comparisons (within sample)

bics
ebics<-exp(-0.5*(bics-min(bics)))
ebics
probs <- round(ebics/sum(ebics), 5)
probs


## Proportion of Deviance  explained, compared to null model
D2s<-c(D.sqr(reg1), D.sqr(reg2), D.sqr(reg3), D.sqr(reg4), D.sqr(reg5),D.sqr(reg6), D.sqr(reg7))
round(D2s, 3)

D2test<- c(Dsquared(reg1), Dsquared(reg2), Dsquared(reg3), Dsquared(reg4),
           Dsquared(reg5), Dsquared(reg6), Dsquared(reg7))

round(D2test,3)

RsqGLM <- c(RsqGLM(reg1), RsqGLM(reg2),RsqGLM(reg3),RsqGLM(reg4),RsqGLM(reg5),
            RsqGLM(reg6),RsqGLM(reg7))

table(RsqGLM)
RsqGLM
## Plots residual diagnostics
models<-list(reg1=reg1, reg2=reg2, reg3=reg3, reg4=reg4, reg5=reg5, reg6=reg6, reg7=reg7)

## randomized quantile residuals
for(i in 1:length(models)){
  r<-models[[i]]
  qr<-qresiduals(r)
  par(mfrow=c(1,2))
  plot(r$fitted, qr, col=as.numeric(dat2$pf_pa)+1, main=names(models)[i], xlab="fitted",
       ylab="randmized quantile residuals")
  abline(h=0)
  qqnorm(qr)
  abline(0,1, col=2)
}


qr5<-qresiduals(reg5)
length(qr)
rf<- fitted.values(reg5)
length(rf)
par(mfrow=c(1,2))
plot(rf, qr5, col=as.numeric(dat2$pf_pa)+1, ylim=c(-4,4), xlab="fitted",
     ylab="Randomized quantile residuals", main = expression(atop('Prevalence '*italic(P.~falciparum), ' Africa')),
     cex.main=1.2, font.main=1)
abline(h=0)
legend(0.11, 4.7, "a", cex = 2.5, xpd=TRUE, bty="n")
qqnorm(qr, main="Normal Q-Q Plot \nAfrica", ylim=c(-4,4), cex.main=1.2, font.main=1)
abline(0,1, col=2)
legend(1, 4.7, "a", cex=2.5, xpd=TRUE, bty="n")

## Visualizations of the best models.

## Plot of the predicted response (STGZ) for AFRICA and ASIA  

labexp <- expression(Probability~of~infection~italic(P.~falciparum))

par(mfrow=c(1,3))
visreg(reg5, "STGZqt1", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (S(T)>0)), ylab= labexp, cex.lab=1.1,
       legend=FALSE)
legend(0.7, 0.15, 
       c("AFRICA", "ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend(0.1, 1, paste('',"a"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)

## Plot of the predicted response (STGZ) at the median log percent of GDP for AFRICA and ASIA.  

visreg(reg7, "lR0G0.gdp", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(log(S(T)>0%*%GDP)), ylab=expression(probability~of~infection~italic(P.~falciparum)), cex.lab=1.1,
       main=expression(Prevalence~italic(Plasmodium~falciparum)), legend=FALSE)
legend(7, 0.15, 
       c("AFRICA", "ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend(7, 1, paste('',"b"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)


## Visualization of the different combinations

my.GDP<-c(min(dat2$lgdpvl), median(dat2$lgdpvl), max(dat2$lgdpvl))
my.GDP

my.lpop<-c(min(dat2$lpop_den), median(dat2$lpop_den), max(dat2$lpop_den))
my.lpop



## STGZ x lGDP
par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "STGZqt1", scale="response", rug=2, overlay=TRUE, 
           cond=list(lpop_den=my.lpop[i], lgdpvl=my.GDP[j]), 
           ylim=c(0,1),line=list(col=c("dodgerblue"), lwd=3),
           fill=list(col=c("#00C1C980")), legend=FALSE,
           xlab=expression(Prob (S(T)>0)), 
           ylab=expression(probability~of~infection~italic(P.~falciparum)))
    legend("topleft",
           c(paste("pop density = ", round(exp(my.lpop[i])), "K", sep=""),
             paste("GDP =", round(exp(my.GDP[j])), " USD", sep ="")),
           bty="n", y.intersp=1)
  }
}


## population density x lGDP

par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "STGZqt1", scale="response", by="AFR", rug=2, overlay=TRUE, 
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

