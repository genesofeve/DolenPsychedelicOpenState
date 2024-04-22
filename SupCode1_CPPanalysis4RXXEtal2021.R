---
title: "Analysis of CPP Behavior for RXX et al 2021"
output: html_notebook
date: "07/06/2021"
editor_options: 
  chunk_output_type: console
---

```{r init,echo=FALSE}
library(tidyverse)
library(plotly)
library(dplyr)
library(tidyr)
library(reshape2)
library(splines)
library(caret)
theme_set(theme_classic())

#preprocess and format data
rd<-read.csv("data/data07112021.csv",stringsAsFactors=FALSE)
rd<-rd[,-1]
colnames(rd)<-paste(rd[1,],rd[2,],rd[3,],rd[4,],c(rep("training",15),rep("test",length(colnames(rd))-15)),sep=".")
rd<-lapply(rd[5:35,],as.numeric)
str(rd)

rdm<-melt(rd)
str(rdm)
rdm<-rdm[-which(is.na(rdm$value)),]
colnames(rdm)[2]<-"variable"
rdm$age<-sapply(rdm$variable,function(x) unlist(strsplit(x,"[.]"))[1])
rdm$age<-as.numeric(gsub("P","",rdm$age))
rdm$genotype<-sapply(rdm$variable,function(x) unlist(strsplit(x,"[.]"))[2])
rdm$genotype[rdm$genotype=="Beta_arrestin"]<-"BetaArr2KO"
rdm$genotype[rdm$genotype=="Beta_arrestin "]<-"BetaArr2KO"
rdm$genotype<-as.factor(rdm$genotype)
rdm$genotype<-relevel(rdm$genotype,ref="WT")
rdm$treatment<-sapply(rdm$variable,function(x) unlist(strsplit(x,"[.]"))[3])
#rdm$treatment[rdm$treatment=="saline"]<-"Saline
rdm$treatment[rdm$treatment==" Psilocybin"]<-"Psilocybin" 
rdm$treatment[rdm$treatment=="ketanserin+LSD"]<-"LSD+Ketanserin" 

rdm$dpt<-sapply(rdm$variable,function(x) unlist(strsplit(x,"[.]"))[4])
rdm$dpt[grep("NA",rdm$dpt)]<-0
rdm$dpt<-as.numeric(rdm$dpt)
rdm$experiment<-sapply(rdm$variable,function(x) unlist(strsplit(x,"[.]"))[5])
str(rdm)

#plot data
fullP<-ggplot(rdm, aes(age, value,color=treatment, linetype=genotype)) + 
  geom_point() +
  stat_smooth(method = lm,formula = y~ns(x,knots=c(42,98))) + scale_shape_manual(values=c(15,19)) + 
  labs(x = "Age", y = "Social preference score \n (normalized)") +
  #scale_color_manual(values=match_col) +
  theme(legend.position = "right")  

fullP


```


```{r plot_params}
#match_col<-c("mediumvioletred","pink4","darkgoldenrod2","tomato4","rosybrown","black","slateblue4","grey30","darkorange","navy","royalpurple","teal")

cmy = function(c, m, y, alpha, maxColorValue=1){
  if(maxColorValue != 1) { c <- c/maxColorValue; m <- m/maxColorValue; y <- y/maxColorValue }
  c <- 1-c; m <- 1-m; y <- 1-y
  hex <- function(v) substring(rgb(v,0,0),2,3)
  if(!missing(alpha)) alpha <- hex(alpha) else alpha <- ''
  paste0('#',hex(c), hex(m), hex(y), alpha)
}

LSD<-c(8,46,100)
Psilocybin<-c(94, 90, 0)
MDMA<-c(40, 100, 100)
Ketamine<-c(22, 88, 0)
Ibogaine<-c(57, 0, 99)
None<-c(32, 27, 26) 
Saline<-c(32, 27, 26) 

cmyFUN<-function(x){cmy(x[1],x[2],x[3],1,max(x))}

```


```{r Original_Nature_curve}

originalOnly<-ggplot(rdm[rdm$treatment=="None"&rdm$genotype=="WT",], aes(age, value)) + 
  stat_smooth(method = lm,formula = y~ns(x,knots=c(35,98))) + scale_shape_manual(values=c(15,19)) + 
  scale_color_manual(values=cmyFUN(None)) +
  labs(x = "Age", y = "Social preference score \n (normalized)") +
  theme(legend.position = "right")  

originalOnly

pdf("original_natureCurve+p126.pdf")
originalOnly
dev.off()
```


```{r subset_to_basic_conditions}

kIndx<-grep("Ketanserin",rdm$treatment)
bIndx<-which(rdm$genotype=="BetaArr2KO")

rdmW<-rdm[-c(kIndx,bIndx),]
sort(unique(rdmW$treatment))
match_col<-sapply(list(Ibogaine,Ketamine,LSD,MDMA,None,Psilocybin,Saline,Saline), cmyFUN)


standardP<-ggplot(rdmW, aes(age, value,color=treatment)) + 
  stat_smooth(method = lm,formula = y~ns(x,knots=c(35,98))) + scale_shape_manual(values=c(15,19)) + 
  labs(x = "Age", y = "Social preference score \n (normalized)") +
  scale_color_manual(values=match_col) +
  theme(legend.position = "right")  

standardP

```

```{r splineFit2Original}

#subset to original controls only 
rdmc<-rdm[which(rdm$experiment=="training"),]

# Split the data into training and test set randomly 
#set.seed(123)
training.samples <- rdmc$variable %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- rdmc[training.samples, ]
test.data <- rdmc[-c(training.samples), ]

#fit beta-spline model with knots at 42, 98
bmodel <- lm(value~0+bs(train.data$age,Boundary.knots=c(35,98)),data = train.data)
summary(bmodel)
#Residual standard error: 0.7561 on 276 degrees of freedom
#Multiple R-squared:  0.5554,  Adjusted R-squared:  0.5506 
#F-statistic: 114.9 on 3 and 276 DF,  p-value: < 2.2e-16


#fit ns-spline model with knots at 42, 98
train.data$splines <- ns(train.data$age,knots=c(35,98))
model <- lm(value~splines,data = train.data)
summary(model)
#Residual standard error: 0.1957 on 275 degrees of freedom
#Multiple R-squared:  0.1053,  Adjusted R-squared:  0.09556 
#F-statistic: 10.79 on 3 and 275 DF,  p-value: 1.003e-06

# check residuals 
plot (train.data$age ,bmodel$resid,xlab = "Age" , ylab = " Residual " )
plot (model$fitted.values ,bmodel$resid,xlab = "Fitted Value" , ylab = " Residual " )

pdf("plots/residualsVsAge.pdf")
plot (train.data$age ,model$resid,xlab = "Age" , ylab = " Residual " )
abline (h=0)
dev.off()

pdf("plots/residualsVsFit.pdf")
plot (model$fitted.values ,model$resid,xlab = "Fitted Value" , ylab = " Residual " )
abline (h=0)
dev.off()

#add model to df 
train.data<-cbind(train.data,predict(model, interval = 'confidence'))

# Model performance
data.frame(
  RMSE = RMSE(train.data$fit, train.data$value),
  R2 = R2(train.data$fit, train.data$value)
)
#       RMSE        R2
# 0.1942869 0.1053173

pdf("plots/control_curve.pdf")
  ggplot(train.data) + geom_point( aes(age, value),size=.5) +
  geom_line(aes(age, fit),size=3) +
  geom_ribbon(aes(x=age, ymin=lwr,ymax=upr), alpha=0.3) +
  labs(x = "Age", y = "Social preference score \n (normalized)") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme(legend.position = "right") +  geom_hline(yintercept = 1,linetype="dashed") 
dev.off()
```

```{r testFit2newControls_Sup2a}
# use new controls as test data
testing.samples<-intersect(grep("saline",rdm$variable),grep("WT",rdm$variable))
testing.samples<-c(testing.samples,intersect(grep("Saline",rdm$variable),grep("WT",rdm$variable)))
test.data <- rdm[testing.samples, ]

#update rdmc to include new controls 
rdmc<-rbind(rdmc,test.data)

# Make predictions
test.data$splines<-ns(test.data$age,knots=c(35,98))
predictions <- model %>% predict(test.data, interval = 'confidence')
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$value),
  R2 = R2(predictions, test.data$value)
)
#       RMSE        R2
#  0.190984 0.01294398

# use new controls as test data
rdmc<-cbind(rdmc, rbind(predict(model, interval = 'confidence'),predictions))
str(rdmc)

pdf("plots/controls.pdf")
  ggplot(rdmc, aes(age, value,color=treatment)) + #geom_point()+
  stat_smooth(method = lm,formula = y~ns(x,knots=c(35,98))) + scale_shape_manual(values=c(15,19)) + 
  labs(x = "Age", y = "Social preference score \n (normalized)") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme(legend.position = "right")  + geom_hline(yintercept = 1,linetype="dashed")
dev.off()


#treatment effects for fit spline model with knots at 35, 98
rdmW$splines <- ns(rdmW$age,knots=c(35,98))
rdmW$treatment<-relevel(as.factor(rdmW$treatment),"None")
model <- lm(value~treatment+splines,data = rdmW)
summary(model)

```

```{r Post_treatment_lenght_of_effect}
str(rdmW)
table(rdmW$dpt)
rdmW$dpt<-as.factor(rdmW$dpt)
rdmW$dpt<-relevel(rdmW$dpt,ref=0)

model <- lm(value~treatment:dpt,data = rdmW)
summary(model)

```


```{r 5HTR2A_Ketaserin_S3}

rdmK<-rdm[-bIndx,]
rdmK$splines <- ns(rdmK$age,knots=c(35,98))
rdmK$treatment<-relevel(as.factor(rdmK$treatment),"None")
model <- lm(value~treatment+splines,data = rdmK)
summary(model)

```

```{r BArrestinKO_KOcontrol_S4}

rdmB<-rdm[-kIndx,]
rdmB$splines <- ns(rdmB$age,knots=c(35,98))
rdmB$treatment<-relevel(as.factor(rdmB$treatment),"None")
rdmB$genotype<-as.factor(rdmB$genotype)
rdmB$genotype<-relevel(rdmB$genotype,ref="WT")
rdmB$age<-as.factor(rdmB$age)

model <- lm(value~genotype:age+splines,data = rdmB[rdmB$treatment=="None",])
summary(model)

```


```{r BArrestinKO_wtreatment_S5}

rdmB$dpt<-as.factor(rdmB$dpt)
rdmB$dpt<-relevel(rdmB$dpt,ref=0)
model <- lm(value~genotype:treatment:dpt+splines,data = rdmB)
summary(model)


```
