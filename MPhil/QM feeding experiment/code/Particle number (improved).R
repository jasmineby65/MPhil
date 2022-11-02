rm(list = ls())

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
require(interactions)

faeces <- as.data.frame(read.csv("../../data/raw/Chap3/Particle number(faeces).csv", header = TRUE, stringsAsFactors = F))
pseudofaeces <- as.data.frame(read.csv("../../data/raw/Chap3/Particle number(pseudofaeces).csv", header = TRUE, stringsAsFactors = F))

head(faeces)
faeces <- faeces[,-c(8,10)]
nrow(faeces)

head(pseudofaeces)
pseudofaeces <- pseudofaeces[, -c(8,9)]
nrow(pseudofaeces)

data <- rbind(faeces, pseudofaeces)
head(data)

data$situ[data$density == "1"] <- "singleton"
data$situ[data$density == "32" & data$Position == "top"] <- "32t"
data$situ[data$density == "32" & data$Position == "bottom"] <- "32b"
data$situ[data$density == "64" & data$Position == "top"] <- "64t"
data$situ[data$density == "64" & data$Position == "bottom"] <- "64b"

data$density<-as.factor(data$density)
data$situ <- as.factor(data$situ)
data$Position <- as.factor(data$Position)

# data$tissue <- exp(-8.91077)*data$size^2.00359
# data$no <- data$ogno/data$tissue
# data$no <- round(data$no, digits = 0)

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
length(x)/nrow(faeces)


for(i in x){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$Bothzero[j] = 1
    }
  }
}

head(data)

faeces <- subset(data, type == "faeces")
pseudofaeces <- subset(data, type == "pseudofaeces")


## trend of 0 ##
zero1 <- glm(Bothzero ~ density, data=faeces, family=binomial)
summary(zero1)
anova(zero1, test="Chi")

faeces1 <- subset(faeces, density != "1")
faeces1

zero2 <- glm(Bothzero ~ density*Position, data=faeces1, family=binomial)
summary(zero2)
anova(zero2, test="Chi")

zero3 <- glm(Bothzero ~ density + Position, data=faeces1, family=binomial)
summary(zero3)
anova(zero3, test="Chi")


## Faeces ##
faeces0 <- subset(faeces, Bothzero == "0")
ggdensity(faeces0$ogno)
ggdensity(faeces0$no)

model1 <- glm.nb(ogno~density,data=faeces0)
par(mfrow = c(2,2))
plot(model1)
summary(model1)
anova(model1, test = "Chi")

library(multcomp)
summary(glht(model1, mcp(density="Tukey")))

# model2 <- glm.nb(ogno~density+size,data=faeces0)
# plot(model2)
# summary(model2)
# anova(model2, test = "Chi")
# 
# model2 <- glm.nb(no~situ,data=faeces0)
# plot(model2)
# summary(model2)
# 
# anova(model2, test = "Chi")
# 
# interact_plot(model2, pred = size, modx = situ, plot.points = T)

faeces01 <- subset(faeces0, density != "1")
head(faeces01)

model3 <- glm.nb(ogno ~ density*Position, data = faeces01)
summary(model3)
anova(model3)

model4 <- update(model3, ~. -density:Position)
summary(model4)
anova(model3, model4, test = "Chi")
anova(model4)

# model5 <- update(model4, ~. -Position:size)
# summary(model5)
# plot(model5)
# anova(model4, model5, test = "Chi")
# anova(model5)
# 
# interact_plot(model5, pred = size, modx = density, plot.points = T)
# unique(faeces$ogno[faeces$density == "32"])


#figure
data1<-ddply(faeces, c("density","Position"), summarise,
             N    = length(ogno),
             mean = mean(ogno),
             sd   = sd(ogno),
             se   = sd / sqrt(N))
data1
fig1<-ggplot(data=faeces0, aes(y=ogno, x=density,fill=Position)) + 
  #geom_errorbar(data=data1, mapping=aes(ymin=0, ymax=mean+sd),position = position_dodge(0.3),width = 0.2) + 
  geom_boxplot(outlier.shape = NA)+
  #geom_point(position = position_dodge(0.3),size=3,aes(shape=Position))+
  theme_classic() +
  scale_fill_grey(start=0.5,end=1)+
  xlab("Density (number of mussels)")+
  ylab("Number of particle") +
  ggtitle("(a) Faeces") 
fig1



### pseudofaeces ###
# Removing 0 #
pseudofaeces0 <- subset(pseudofaeces, Bothzero == "0")

ggdensity(pseudofaeces0$ogno)
# pseudofaeces0$no <- pseudofaeces0$ogno/pseudofaeces0$size
# pseudofaeces0$density<-as.factor(pseudofaeces0$density)

model1 <- glm.nb(ogno~density,data=pseudofaeces0)
par(mfrow = c(2,2))
plot(model1)
summary(model1)
anova(model1, test = "Chi")

pseudofaeces01 <- subset(pseudofaeces0, density != "1")

model2 <- glm.nb(ogno~density*Position,data=pseudofaeces01)
par(mfrow = c(2,2))
plot(model2)
summary(model2)
anova(model2, test = "Chi")

model3 <- glm.nb(ogno~density+Position,data=pseudofaeces01)
par(mfrow = c(2,2))
plot(model3)
summary(model3)
anova(model3, test = "Chi")

#figure
data3<-ddply(pseudofaeces, c("density","Position"), summarise,
             N    = length(no),
             mean = mean(no),
             sd   = sd(no),
             se   = sd / sqrt(N))

fig3<-ggplot(data=pseudofaeces0, aes(y=ogno, x=density,fill=Position)) + 
  geom_boxplot(outlier.shape = NA)+
  theme_classic() +
  scale_fill_grey(start=0.5,end=1)+
  xlab("Density (number of mussels)")+
  ylab("Number of particle/ml") +
  ggtitle("(b) Pseudofaeces") + ylim(0,9)
fig3

pdf("../../results/mono/Chap3/particle_no_mono.pdf", width = 6, height = 11)
grid.arrange(fig1,fig3)
dev.off()
