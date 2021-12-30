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

##Combinations

dat2$l.prev<-log(dat2$pf_pa*dat2$pf_pr+1)

##scaled
##dat2$lpop_den01<-dat2$lpop_den/max(dat2$lpop_den)
##dat2$lGDP01<-dat2$lGDP/max(dat2$lGDP)

## we consider a few combinations of R0 measures and population
summary(dat2$Roq1)
summary(dat2$pop_den)

dat2$R0G0.lpop<-dat2$R0G0qt1*log(dat2$pop_den)
dat2$R0G0.lgdp<- dat2$R0G0qt1*log(dat2$gdpvl)
dat2$lR0.pop<- log((dat2$Roq1*dat2$pop_den)+1)
dat2$lR0G0.pop<-log((dat2$R0G0qt1*dat2$pop_den)+1)
dat2$lR0.gdp<- log((dat2$Roq1*dat2$gdpvl)+1)
dat2$lR0G0.gdp<- log((dat2$R0G0qt1*dat2$gdpvl)+1)

summary(dat2$R0G0.lpop)
summary(dat2$R0G0.lgdp)
summary(dat2$lR0.pop)
summary(dat2$lR0G0.pop)
summary(dat2$lR0.gdp)
summary(dat2$lR0G0.gdp)


#d$R0.GR0.alt<-d$R0.GR0*log(d$population)
## Preliminary data visualization

#First we visually explore our covariates. Notice there is quite a bit of correlation (colinearity) between our environmental (T.lag, R0) and socio-economic variables, and between the socioeconomic variables, especally tourism GDP% and log population. This isn't causal, but this means we have to take care when interpretting the coefficients in the regressions
dat2$temp <- dat2$qtmean1
hist(dat2$temp)

names(dat2)
dim(dat2)




vary<- c(2,6,9,10,11,12,17,19)
pairs(dat2[,vary])


round(cor(dat2[,c(2,6,9,10,11,12,17,19)]), digits=2)
names(dat2)
## Models for the Presence/Absence Data


## 1 lpop_den only
reg1<-glm(pf_pa ~ lpop_den, family=binomial, data=dat2)
summary(reg1)

## 2 lGDP only
reg2<-glm(pf_pa ~ lgdpvl, family=binomial, data=dat2)
summary(reg2)

## 3 ROGO only
reg3<-glm(pf_pa ~ R0G0qt1, family=binomial, data=dat2)
summary(reg3)

## 4 socio factors only
reg4<-glm(pf_pa ~ (lgdpvl +lpop_den), family=binomial, data=dat2)
summary(reg4)

## 5 ROGO + lGDP + lpop_den
reg5<-glm(pf_pa ~ (R0G0qt1+lgdpvl+lpop_den), family=binomial, data=dat2)
summary(reg5)

## 6 log(ROGO*lpop_den)
reg6<-glm(pf_pa ~ lR0G0.pop, family=binomial, data=dat2)
summary(reg6)

## 7 log(ROGO*gdp)
reg7<-glm(pf_pa ~ lR0G0.gdp, family=binomial, data=dat2)
summary(reg7)


bics <- c(BIC(reg1), BIC(reg2), BIC(reg3), BIC(reg4), BIC(reg5), BIC(reg6), BIC(reg7)) 


#### Model Comparisons (within sample)

bics
ebics<-exp(-0.5*(bics-min(bics)))
ebics
probs <- round(ebics/sum(ebics), 5)
probs

summ(reg5, confint=TRUE, digits=3)

summary(reg5)
## Quantitative methods in Ecology and Evolution , Biology 548b
## Comparison of estimators for SeriesB Series
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
legend(-6, 4.7, "b", cex=2.5, xpd=TRUE, bty="n")

#Reg5.sub is the top model, by BIC, for the full data set. We re-fit with all the predictors scaled,
#so we can see the relative size of the parameters for important parameters. 
#Note that due to multi-colinearity in the predictors,
#we can't reliably interpret these parameters as relating to the relative influence of various predictors (or even as the marginal impact of a predictor).

summary(glm(pf_pa ~ (R0G0qt1+lGDP01+lpop_den01),family=binomial, data=dat2))


### Visualizations of the model fit.

##Here are some visualizations of the model fit. First I plot the predicted response at the median log percent of GDP,  across a subset of the values of proportion of GDP in toursim and across log population from the data set. The two sets of lines indicate the fit for Dengue (DEN: 1) and for the other 2 diseases (DEN: 0)


##png(file="PA_reg1_cond_resp_GDP7274.png", height=900, width=900, pointsize=18)
labexp <- expression(Probability~of~infection~italic(P.~falciparum))

par(mfrow=c(1,3))
visreg(reg5, "R0G0qt1", scale="response", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (S(T)>0)), ylab= labexp, cex.lab=1.1,
       legend=FALSE)
legend(-0.18, 1.05, paste('',"a"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)


visreg(reg5, "lgdpvl", scale="response", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab= "log(GDP)", ylab= labexp, cex.lab=1.1,
       legend=FALSE)
legend(4.8, 1.05, paste('',"b"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)

visreg(reg5, "lpop_den", scale="response", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab="log(population density)", ylab= labexp, cex.lab=1.1,
       legend=FALSE)
legend(-0.15, 1.05, paste('',"c"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)



visreg(reg7, "lR0G0.gdp", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(log(S(T)>0%*%GDP)), ylab=expression(probability~of~infection~italic(P.~falciparum)), cex.lab=1.1,
       main=expression(Prevalence~italic(Plasmodium~falciparum)), legend=FALSE)
legend(7, 0.15, 
       c("AFRICA", "ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend(-1.7, 1.04, paste('',"b"), text.font = 3, cex = 2.5,
       bty="n", y.intersp=0.9)


#"#FF4E3780","#00C1C980"
#line=list(col=c("red","dodgerblue"), lwd=3),
#fill=list(col=c("red","dodgerblue")), 
#scale="response",


summary(dat2$lgdpvl)
summary(dat2$lpop_den)

my.GDP<-c(min(dat2$lgdpvl), median(dat2$lgdpvl), max(dat2$lgdpvl))
my.GDP

my.lpop<-c(min(dat2$lpop_den)-1, median((dat2$lpop_den)-5), max((dat2$lpop_den)-5))
my.lpop



##Combinations
## temperature x lGDP
par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "R0G0qt1", scale="response", rug=2, overlay=TRUE, 
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


##population density vs GDP

par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "R0G0qt1", scale="response", by="AFR", rug=2, overlay=TRUE, 
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



par(mfrow=c(3,3))
for (i in 1:3){
  for(j in 1:3){
    visreg(reg5, "R0G0qt1", scale="response", by="AFR", rug=2, overlay=TRUE, 
           cond=list(lpop_den=my.lpop[i], lGDP=my.GDP[j]), 
           ylim=c(0,1),line=list(col=c("red","dodgerblue"), lwd=3),
           fill=list(col=c("#FF4E3780","#37b6ff")), legend=FALSE,
           xlab=expression(Prob (S(T)>0)), ylab="probability prevalence P. falciparum")
    legend("topleft",
           c(paste("pop density = ", round(10^(my.lpop[i])), "K", sep=""),
             paste("GDP =", round(10^my.GDP[j]), " USD", sep ="")),
           bty="n", y.intersp=1)
  }
}



##SECOND PART HURDLE MODEL
## POSITIVE VALUES ONLY

dd<-subset(dat2, pa!=0)
names(dd)
dim(dd)
summary(dd)
hist(dd$pf_pr)


##MODELS

## 1 lpop_den only
mp1<-glm(pf_pr ~ lpop_den*AFR - AFR, family="Gamma", data=dd)
summary(mp1)

## 2 lGDP only
mp2<-glm(pf_pr ~ lGDP*AFR - AFR, family="Gamma", data=dd)
summary(mp2)

## 3 socio factors only
mp3<-glm(pf_pr ~ (lGDP +lpop_den)*AFR - AFR, family="Gamma", data=dd)
summary(mp3)

## 4 R0GO
mp4<-glm(pf_pr ~ R0G0qt1*AFR - AFR, family="Gamma", data=dd)
summary(mp4)

##5 R0GO + lpop_den
mp5<-glm(pf_pr ~ (R0G0qt1+lpop_den+lGDP)*AFR - AFR, family="Gamma", data=dd)
summary(mp5)

##6 R0GO * lgdp
mp6<-glm(pf_pr ~ (lR0G0.gdp)*AFR - AFR, family="Gamma", data=dd)
summary(mp6)


bics <- c(BIC(mp1), BIC(mp2), BIC(mp3),BIC(mp4), BIC(mp5), BIC(mp6)) 

#### Model Comparisons (within sample)

bics
ebics<-exp(-0.5*(bics-min(bics)))
probs <- round(ebics/sum(ebics), 5)
probs

## Proportion of Deviance  explained, compared to null model
D2s<-c(D.sqr(mp1), D.sqr(mp2),  D.sqr(mp3),D.sqr(mp4), D.sqr(mp5),D.sqr(mp6))

round(D2s, 3)


## Plots residual diagnostics
models<-list(mp1=mp1, mp2=mp2, mp3=mp3, mp4=mp4, mp5=mp5, mp6=mp6)

## randomized quantile residuals
par(mfrow=c(1,2))

for(i in 1:length(models)){
  r<-models[[i]]
  qr<-qresiduals(r)
  par(mfrow=c(1,2))
  plot(r$fitted, qr, col=as.numeric(dd$pf_pr)+1, main=names(models)[i], xlab="fitted")
  abline(h=0)
  qqnorm(qr)
  abline(0,1, col=2)
}



qr<-qresiduals(mp5)
par(mfrow=c(1,2))
plot(r$fitted, qr, col=as.numeric(dd$pf_pr)+1, xlab="fitted", 
     ylab = "randomized quantile residuals")
abline(h=0)
qqnorm(qr)
abline(0,1, col=2)


#Reg1.sub is the top model, by BIC, for the full data set. We re-fit with all the predictors scaled,
#so we can see the relative size of the parameters for important parameters. 
#Note that due to multi-colinearity in the predictors,
#we can't reliably interpret these parameters as relating to the relative influence of various predictors (or even as the marginal impact of a predictor).

summary(glm(pf_pr ~ R0G0qt1+(lGDP01+lpop_den01),family="Gamma", data=dd))


### Visualizations of the model fit.

##Here are some visualizations of the model fit. First I plot the predicted response at the median log percent of GDP,  across a subset of the values of proportion of GDP in toursim and across log population from the data set. The two sets of lines indicate the fit for Dengue (DEN: 1) and for the other 2 diseases (DEN: 0)


##png(file="PA_reg1_cond_resp_GDP7274.png", height=900, width=900, pointsize=18)

par(mfrow=c(1,1))
visreg(mp5, "R0G0qt1", scale="response", rug=2, by="AFR", overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (R0>0)), ylab="probability prevalence P. vivax", 
       main="Prevalence Plasmodium vivax", legend=FALSE)
legend("topleft", 
       c("AFRICA","ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue","red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend("topright", paste('',"B"), text.font = 3, cex = 1.8,
       bty="n", y.intersp=0.9)


par(mfrow=c(1,1))
visreg(mp5, "lGDP", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (lGDP)), ylab="probability prevalence P. vivax", 
       main="Prevalence Plasmodium vivax", legend=FALSE)
legend("topleft", 
       c("AFRICA","ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue","red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend("topright", paste('',"B"), text.font = 3, cex = 1.8,
       bty="n", y.intersp=0.9)


par(mfrow=c(1,1))
visreg(mp5, "lpop_den", scale="response", by="AFR", rug=2, overlay=TRUE, ylim=c(0,1),
       xlab=expression(Prob (R0>0)), ylab="probability prevalence P. vivax", 
       main="Prevalence Plasmodium vivax", legend=FALSE)
legend("topleft", 
       c("AFRICA","ASIA"), cex = 1.2, text.font =3,
       col=c("dodgerblue", "red"), lwd=c(3,3),
       bty="n", y.intersp=1)
legend("topright", paste('',"B"), text.font = 3, cex = 1.8,
       bty="n", y.intersp=0.9)



