---
title: "Manova"
output: html_notebook
date: "07/06/2021"
editor_options: 
  chunk_output_type: console
---

```{r init,echo=FALSE
library(MASS)
library(mvnormtest)
library(reshape)
library(caret)
library(ggplot2)

dat<-read.csv("data/ephys.csv")
str(dat)

dat$Age[dat$Age==1]<-"P42"
dat$Age[dat$Age==2]<-"P98"
dat$Structure[dat$Structure==1]<-"Nac"
dat$Structure[dat$Structure==2]<-"mPFC"
dat$treatment[dat$treatment==1]<-"sal"
dat$treatment[dat$treatment==2]<-"MDMA"
dat$treatment[dat$treatment==3]<-"LSD"
dat$treatment[dat$treatment==4]<-"psilo"
dat$treatment[dat$treatment==5]<-"keta"
dat$treatment[dat$treatment==6]<-"Ibo"

```

```{r manova}

sIndx<-dat$treatment=="sal"
paste0(dat$treatment)


datM<-melt(dat, id=c('Age','Structure','session','treatment'),measured=c('frequency','amplitute'))
str(datM)

datM<-melt(dat, id=c('treatment'),measured=c('frequency','amplitute'))
str(datM)

dat$treatment<-as.factor(dat$treatment)
dat$treatment<-relevel(dat$treatment,ref="sal")
dat$Structure<-as.factor(dat$Structure)
dat$Structure<-relevel(dat$Structure,ref="Nac")
dat$session<-as.factor(dat$session)

mlmFull<-lm(cbind(dat$frequency,dat$amplitute)~treatment:Age*Structure,data=dat)
mlm<-lm(cbind(dat$frequency,dat$amplitute)~Age*Structure,data=dat)

anova(mlmFull,mlm)

```
