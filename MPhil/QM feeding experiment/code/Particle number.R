library(lme4)
library(Rcpp)
library(ggplot2)
library(ggpubr)
library(plyr)
library(pscl)
library(MASS)
library(glmmTMB)
library(gridExtra)
library(lsmeans)
library(lmtest)
library(multcomp)

data <- read.csv(file.choose(), header = TRUE)
head(data)
data$density<-as.factor(data$density)

#Trend of 0#
#Finding individuals with no capsule in both faeces and pseudofaeces
data$Bothzero <- 0
head(data)

nrow(data)
unique(data$type)

z <- vector()
for(i in 1:nrow(data)){
  if(data$zero[i] == "1" && data$type[i] == "faeces"){
  z <- c(z, data$ID[i])  
  }
}
z

y <- vector()
for(i in 1:nrow(data)){
  if(data$zero[i] == "1" && data$type[i] == "pseudofaeces"){
    y <- c(y, data$ID[i])  
  }
}
y

x <- intersect(z,y)
x

for(i in x){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$Bothzero[j] = 1
    }
  }
}

head(data)


zero1 <- glm(zero ~ density+Position, data=faeces, family=binomial)
summary(zero1)

zero2 <- glm(zero ~ Position, data=faeces, family=binomial)
summary(zero2)
anova(zero1, zero2, test="Chi")

zero3 <- glm(zero ~ density, data=faeces, family=binomial)
summary(zero3)
anova(zero1, zero3, test="Chi")


###Faeces
#trend of 0
faeces <- read.csv(file.choose(), header = TRUE)
head(faeces)
faeces$density<-as.factor(faeces$density)

zero1 <- glm(zero ~ density+Position, data=faeces, family=binomial)
summary(zero1)

zero2 <- glm(zero ~ Position, data=faeces, family=binomial)
summary(zero2)
anova(zero1, zero2, test="Chi")

zero3 <- glm(zero ~ density, data=faeces, family=binomial)
summary(zero3)
anova(zero1, zero3, test="Chi")


# Removing 0


#zer-inflated model
#without singletons for position analysi
faeces1 <- read.csv(file.choose(), header = TRUE)
faeces1$density<-as.factor(faeces1$density)
model3<- zeroinfl(ogno ~ density+Position+size, dist="negbin", data = faeces1)
summary(model3)
model4<- zeroinfl(ogno ~ size*Position, dist="negbin", data = faeces1)
summary(model4)
model4$df.residual
drop1(model4, test="Chisq")

#with 1 but no position
faeces$density<-as.factor(faeces$density)
model5<- zeroinfl(ogno ~ size*density, dist="negbin", data = faeces)
summary(model5)
model5$df.residual
model6<- zeroinfl(ogno ~ density+ size, dist="negbin", data = faeces)
summary(model6)
model6$df.residual
drop1(model6,test="Chisq")
anova(model4,model5)
faeces64 <- read.csv(file.choose(), header = TRUE)
faeces64$density<-as.factor(faeces64$density)
model7<- zeroinfl(ogno ~ density, dist="negbin", data = faeces64)
summary(model7)

#no 0 data
faeces0 <- read.csv(file.choose(), header = TRUE)
ggdensity(faeces0$ogno)
faeces0$density<-as.factor(faeces0$density)

plot(ogno~size,data=faeces0)

model0.3 <- glm.nb(ogno~density+Position,data=faeces0)
summary(model0.3)
faeces0.1 <- within(faeces0, density <- relevel (faeces0$density, ref="32"))
model0.2 <- glm.nb(ogno~Position+density,data=faeces0.1)
summary(model0.2)

faeces01 <- read.csv(file.choose(), header = TRUE)
ggdensity(faeces0$ogno)
faeces01$density<-as.factor(faeces01$density)
model0.4 <- glm.nb(ogno~density+Position,data=faeces01)
summary(model0.4)


#figure
data1<-ddply(faeces, c("density","Position"), summarise,
             N    = length(no),
             mean = mean(no),
             sd   = sd(no),
             se   = sd / sqrt(N))

fig1<-ggplot(data=faeces0, aes(y=no, x=density,fill=Position)) + 
  #geom_errorbar(data=data1, mapping=aes(ymin=0, ymax=mean+sd),position = position_dodge(0.3),width = 0.2) + 
  geom_boxplot(outlier.shape = NA)+
  #geom_point(position = position_dodge(0.3),size=3,aes(shape=Position))+
  theme_classic() +
  scale_fill_grey(start=0.5,end=1)+
  xlab("Density (number of mussels)")+
  ylab("Number of particle/g tissue dry mass") +
  ggtitle("(a) Faeces") 
fig1


###pseudofaeces
#trend of 0
pseudofaeces <- read.csv(file.choose(), header = TRUE)
pseudofaeces$density<-as.factor(pseudofaeces$density)
zero3 <- glm(zero ~Position+density+size, data=pseudofaeces, family=binomial)
summary(zero3)
pseudofaeces1 <- read.csv(file.choose(), header = TRUE)
pseudofaeces1$density<-as.numeric(pseudofaeces1$density)
zero4 <- glm(zero ~Position+size, data=pseudofaeces1, family=binomial)
summary(zero4)

#zer-inflated model
#without singletons for position analysis
pseudofaeces1$density<-as.numeric(pseudofaeces1$density)
model8<- zeroinfl(ogno ~ density+Position+size, dist="negbin", data = pseudofaeces1)
summary(model8)
model8$df.residual

#with 1 but no position
pseudofaeces$density<-as.factor(pseudofaeces$density)
model9<- zeroinfl(ogno ~ size+density, dist="negbin", data = pseudofaeces)
summary(model9)
model9$df.residual

#no 0 data
pseudofaeces0 <- read.csv(file.choose(), header = TRUE)
ggdensity(pseudofaeces0$ogno)
pseudofaeces0$density<-as.factor(pseudofaeces0$density)

model0.5 <- glm.nb(ogno~density+Position+size,data=pseudofaeces0)
summary(model0.5)
pseudofaeces0.1 <- within(pseudofaeces0, density <- relevel (pseudofaeces0$density, ref="32"))
model0.6 <- glm.nb(ogno~Position+density,data=pseudofaeces0.1)
summary(model0.6)

pseudofaeces01 <- read.csv(file.choose(), header = TRUE)
ggdensity(faeces0$ogno)
pseudofaeces01$density<-as.factor(pseudofaeces01$density)
model0.7 <- glm.nb(ogno~density+Position,data=pseudofaeces01)
summary(model0.7)

#figure
data3<-ddply(pseudofaeces, c("density","Position"), summarise,
             N    = length(no),
             mean = mean(no),
             sd   = sd(no),
             se   = sd / sqrt(N))

fig3<-ggplot(data=pseudofaeces0, aes(y=no, x=density,fill=Position)) + 
  geom_boxplot(outlier.shape = NA)+
  theme_classic() +
  scale_fill_grey(start=0.5,end=1)+
  xlab("Density (number of mussels)")+
  ylab("Number of particle/ml/g tissue dry mass") +
  ggtitle("(b) Pseudofaeces") 
fig3
grid.arrange(fig1,fig3)
