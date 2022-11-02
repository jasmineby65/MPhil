###########################
###### Particle size ######
###########################

rm(list = ls())

library(lme4)
library(Rcpp)
library(ggplot2)
library(ggpubr)
library(lsmeans)
library(gridExtra)

### Faeces & pseudofaeces ###
data <- as.data.frame(read.csv("../../data/raw/Chap3/Particle size.csv", header = TRUE, stringsAsFactors = F, na.strings=c("","NA")))
head(data)

data <- data[,-c(5,7)]
data$type[data$density == "control"] <- "control"
data <- subset(data, type != "control")
data$size <- as.integer(data$size)
data$individual[is.na(data$individual)] <- data$batch[is.na(data$individual)]
head(data)

data$situ <- NA
data$situ[data$density == "singleton"] <- "singleton"
data$situ[data$density == 32 & data$position == "top"] <- "32t"
data$situ[data$density == 32 & data$position == "bottom"] <- "32b"
data$situ[data$density == 64 & data$position == "top"] <- "64t"
data$situ[data$density == 64 & data$position == "bottom"] <- "64b"

ggdensity(data$size)

model1<-lmer(size~situ+type+(1|density), data=data)
par(mfrow = c(2, 2))
plot(model1)
summary(model1)
anova(model1, test = "Chi")
TukeyHSD((aov(model1)))

model1<-glmer(size~situ+type +(1|density), data=data, family = "Gamma")
model2<-glmer(size~situ +(1|density), data=data, family = "Gamma")
anova(model1, model2)


### Faeces ###
faeces <- as.data.frame(read.csv("../../data/raw/Chap3/Particle size (faeces).csv", header = TRUE, stringsAsFactors = F, na.strings=c("","NA")))
head(faeces)

faeces <- faeces[,-c(4,5,7,8)]
faeces$size <- as.integer(faeces$size)
faeces$situ <- NA
faeces$situ[faeces$density == "control"] <- "control"
faeces$situ[faeces$density == "singleton"] <- "singleton"
faeces$situ[faeces$density == 32 & faeces$Position == "top"] <- "32t"
faeces$situ[faeces$density == 32 & faeces$Position == "bottom"] <- "32b"
faeces$situ[faeces$density == 64 & faeces$Position == "top"] <- "64t"
faeces$situ[faeces$density == 64 & faeces$Position == "bottom"] <- "64b"
head(faeces)

faeces$density <- as.factor(faeces$density)
faeces$batch <- as.factor(faeces$batch)
faeces$Position <- as.factor(faeces$Position)
str(faeces)

ggdensity(faeces$size)

model4 <- glmer(size ~ situ + (1|density), data = faeces, family = "Gamma")
summary(model4)
anova(model4, test="Chisq")

model5 <- glmer(size ~ (1|density), data = faeces, family = "Gamma")
anova(model4, model5)

library(multcomp)
summary(glht(model4, mcp(situ="Tukey")))


model4 <- glm(size~density,data=faeces,family = "Gamma")
plot(model4)
summary(model4)
anova(model4, test="Chisq")

faeces_layer <- subset(faeces, !is.na(faeces$Position))

model5 <- glm(size~density*Position,data=faeces,family = "Gamma")
anova(model5, test="Chisq")
TukeyHSD(aov(model5))

par(mfrow = c(1, 1))
ggplot(data,aes(y=size, x=position, colour=type))+
  geom_boxplot()+
  scale_fill_grey(start=0.5,end=1)+
  theme_classic() +
  ylab("Particle size (??m)")+
  ylim(0,200)+
  xlab("Aggregate density (number of mussels)")+
  ggtitle("(a)")


faeces$density <- factor(faeces$density, levels = c('control','singleton','32','64'),ordered = TRUE)

fig1<-ggplot(faeces,aes(y=size, x=density, fill=Position))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_grey(start=0.5,end=1)+
  theme_classic() +
  ylab("Particle size (?m)")+
  ylim(0,200)+
  xlab("Density (number of mussels)")+
  ggtitle("(a) Faeces")
fig1


### Pseudofaeces ###
pseudofaeces <- as.data.frame(read.csv("../../data/raw/Chap3/Particle size (pseudofaeces).csv", 
                                       header = TRUE, stringsAsFactors = F, na.strings=c("","NA")))
head(pseudofaeces)

pseudofaeces <- pseudofaeces[,-c(4,5,7,8)]
pseudofaeces$size <- as.integer(pseudofaeces$size)
head(pseudofaeces)

pseudofaeces$situ <- NA
pseudofaeces$situ[pseudofaeces$density == "control"] <- "control"
pseudofaeces$situ[pseudofaeces$density == "singleton"] <- "singleton"
pseudofaeces$situ[pseudofaeces$density == 32 & pseudofaeces$position == "top"] <- "32t"
pseudofaeces$situ[pseudofaeces$density == 32 & pseudofaeces$position == "bottom"] <- "32b"
pseudofaeces$situ[pseudofaeces$density == 64 & pseudofaeces$position == "top"] <- "64t"
pseudofaeces$situ[pseudofaeces$density == 64 & pseudofaeces$position == "bottom"] <- "64b"
str(pseudofaeces)
ggdensity(pseudofaeces$size)

model5 <- glm(size ~ density, data = pseudofaeces, family = "Gamma")
par(mfrow = c(2,2))
plot(model5)
summary(model5)
anova(model5, test="Chisq")
TukeyHSD(aov(model5))

pseudo_layer <- subset(pseudofaeces, !is.na(pseudofaeces$position))
pseudo_layer

model6 <- glm(size ~ density * position, data = pseudo_layer, family = "Gamma")
plot(model6)
summary(model6)
anova(model6, test="Chisq")
TukeyHSD(aov(model6))

model6 <- glm(size ~ position + density, data = pseudofaeces, family = "Gamma")
plot(model6)
summary(model6)
anova(model6, test="Chisq")
TukeyHSD(aov(model6))

pseudofaeces <- read.csv(file.choose(), header = TRUE)
pseudofaeces$density <- factor(pseudofaeces$density, levels = c('control','singleton','32','64'),ordered = TRUE)
fig2<-ggplot(pseudofaeces,aes(y=size, x=density, fill=position))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_grey(start=0.5,end=1)+
  theme_classic() +
  ylab("Particle size (?m)")+
  xlab("Density (number of mussels)")+
  ylim(0,200)+
  ggtitle("(b) Pseudofaeces")
fig2
grid.arrange(fig1,fig2)


### Particle abundance ###
data <- as.data.frame(read.csv("../data/raw/Chap3/Particle size.csv", header = TRUE, stringsAsFactors = F, na.strings=c("","NA")))

unique(data$batch)
a <- sum(data$batch == "control1")
a
b <- sum(data$batch == "control3")
b
((a+b)/2)*1000
