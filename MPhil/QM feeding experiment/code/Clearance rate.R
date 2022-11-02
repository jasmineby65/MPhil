rm(list = ls())

library(gridExtra)
library(ggplot2)
library(lsmeans)
library(ggpubr)

trial1 <- as.data.frame(read.csv("../../data/raw/Chap3/Clearance rate trial1.csv", header = TRUE, stringsAsFactors = F))
head(trial1)
trial1 <- trial1[,-c(3,6,7,8)]

trial2 <- as.data.frame(read.csv("../../data/raw/Chap3/Clearance rate trial2.csv", header = TRUE, stringsAsFactors = F))
trial2
trial2 <- trial2[,-c(3:6,8:10,12:19)]
trial2

# Settlement of algae #
control1 <- subset(trial1, density == "control")
control2 <- subset(trial2, density == "control")
trial1 <- rbind(trial1, control2)
control <- rbind(control1, control2)
control

settlement <- (control$start - control$end)/control$start
settlement
mean(settlement)


# Removal #
removal <- trial1
removal$density[removal$density == "control"] <- 0
removal$reduction <- removal$start - removal$end
ggdensity(removal$reduction)
removal$density <- as.factor(removal$density)

model <- lm(reduction ~ density, data = removal) 
anova(model)
summary(model)
TukeyHSD((aov(model)))
model2 <- segmented(model, seg.Z = ~density)
summary(model2)
davies.test(model2)

ggplot(removal, aes(x=density, y=reduction)) +
  # geom_smooth(method = "glm", method.args = list(family = "Gamma"), se=F, colour="grey") + 
  geom_boxplot(aes(group=density)) +
  xlab("Clump density (number of mussels)") + 
  ylab("Clearance rate (ml/mussel/hr)") + 
  theme_classic()


plot(reduction ~ density, data = removal)

# CR calculation #
trial1 <- subset(trial1, density != "control")
trial1$density <- as.numeric(trial1$density)
trial1$CR <- (1000/(trial1$density*2))*(log(trial1$start/trial1$end))
trial1

ggdensity(trial1$CR)
  
model7 <- glm(CR ~ density, family="Gamma", data=trial1)
par(mfrow = c(2,2))
plot(model7)
summary(model7)
anova(model7, test = "Chi")

pdf("../../results/mono/CR.pdf", width = 8, height = 5)
ggplot(trial1, aes(x=density, y=CR)) +
  geom_smooth(method = "glm", method.args = list(family = "Gamma"), se=F, colour="grey") + 
  geom_boxplot(aes(group=density)) +
  xlab("Clump density (number of mussels)") + 
  ylab("Clearance rate (ml/mussel/hr)") + 
  ylim(0, 700) +
  theme_classic()
dev.off()


library(plyr)
ddply(trial1, ~density, dplyr::summarise,
      N    = length(CR),
      mean = mean(CR),
      sd   = sd(CR),
      se   = sd / sqrt(N))

