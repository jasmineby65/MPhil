############################################
### Body condition vs clump size (Chap3) ###
############################################

rm(list = ls())

require(gridExtra)
require(ggplot2)
require(smatr)
require(lme4)
require(Rcpp)
require(lmerTest)
require(plyr)
require(effects)


qmall <- as.data.frame(read.csv("../../data/raw/Chap3/All mussels.csv", header = TRUE, stringsAsFactors = F))
head(qmall)

qmall <- qmall[,-c(4, 5, 6, 9, 10, 11, 12, 14, 15, 16, 17, 19)]
head(qmall)

colnames(qmall)[6] <- "tissue"
colnames(qmall)[5] <- "shell"
qmall$clump <- as.factor(qmall$clump)
qmall$body <- qmall$tissue/qmall$shell
head(qmall)


### Effect of density on size ### 
model1 <- lmer(size ~ density + (1|clump), qmall)
summary(model1)


### SMA calculation ###
## Tissue ##
# ln-transformation
qmall$lntissue <- log(qmall$tissue) 
qmall$lnsize <- log(qmall$size) 

# Fitting SMA regression 
qm_sma1 <- sma(lntissue ~ lnsize, data = qmall) 
qm_sma1


## Shell ## 
# ln-transformation
qmall$lnshell <- log(qmall$shell) 

# Fitting SMA regression 
qm_sma2 <- sma(lnshell ~ lnsize, data = qmall) 
qm_sma2



### clump density vs body condition ###
# All #
model2 <- lmer(body ~ density + size + (1|clump), data=qmall)
summary(model2)

# 26-30mm #
qm <- subset(qmall, size <= 30 & size >= 26)
head(qm)

# Density of living mussels
model5 <- lmer(body ~ density + size + (1|clump), data=qm)
summary(model5)

model6 <- lmer(body ~ density + (1|clump), data=qm)
summary(model6)


# Plotting #
plot <- as.data.frame(effects::effect(term= "density", mod= model6))
head(plot) 


pdf("../../results/mono/Chap3/Condition_density_mono.pdf", width = 8, height = 5)

ggplot() + 
  geom_boxplot(data=qm, aes(x=density,y=body,group=density),outlier.shape = NA, colour="grey") +
  geom_line(data = plot, aes(x = density, y = fit)) +
  geom_ribbon(data = plot, aes(x = density, ymin = lower, ymax = upper), alpha = 0.3) +
  xlab("Clump size (number of mussels per clump)") + 
  ylab("Tissue mass/shell mass (g)") + 
  theme_classic()

dev.off()


### Length-dry weight relationship ### 
eq1 <- lm(log(tissue)~log(size),data=qmall)
summary(eq1)
plot(tissue ~ size , data = qmall)
