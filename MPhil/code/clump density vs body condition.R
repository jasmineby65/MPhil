############################################
####### Body condition vs clump size #######
############################################

rm(list = ls())

require(gridExtra)
require(ggplot2)
require(lsmeans)
require(segmented)
require(interactions)
require(jtools)
require(grid)
require(gridtext)


#######################
#### Quagga mussel ####
#######################

qm <- as.data.frame(read.csv("../data/standardized_QM.csv", header = TRUE, stringsAsFactors = F))
head(qm)

qm$clump.size <- as.numeric(qm$clump.size)

sum(qm$clump.size == 1)

# Fitting segmented model to linear model
qm_model1 <- lm(body ~ clump.size*size, data = qm)
summary(qm_model1)

qm_model2 <- segmented(qm_model1, seg.Z = ~clump.size)
summary(qm_model2)
anova(qm_model2)
davies.test(qm_model1, seg.Z=~clump.size)

slope(qm_model2)
intercept(qm_model2)
confint(qm_model2)

plot(qm_model2)

# Trend of shell length with clump size
qm_model <- lm(size ~ clump.size, data=qm)
summary(qm_model)
anova(qm_model)
qm_model_seg <- segmented(qm_model, seg.Z = ~clump.size)
summary(qm_model_seg)
davies.test(qm_model_seg)


## Plotting ## 
# Interaction effect - colour #
pdf("../results/colour/qm_interaction_colour.pdf", width = 5, height = 5)
interact_plot(qm_model2, pred = clump.size, modx = size,
              x.label = "Aggregate size", y.label = "tissue mass/shell mass (g)", legend.main = "Shell length (mm)") +
  theme(legend.position = c(0.8, 0.18), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) 
dev.off()


# Interaction effect - mono #
pdf("../results/mono/qm_interaction_mono.pdf", width = 5, height = 5)
interact_plot(qm_model2, pred = clump.size, modx = size, colors = "Greys",
              x.label = "Aggregate size", y.label = "tissue mass/shell mass (g)", legend.main = "Shell length (mm)") +
  theme(legend.position = c(0.8, 0.18), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) 
dev.off()


# Body condition vs aggregate density # 
# Predicted values for constant shell length 

fitted <- data.frame(clump.size = seq(1, 310, 1))
fitted$size <- rep(mean(qm$size), nrow(fitted))
head(fitted)

fitted1 <- predict(qm_model2, newdata = fitted, se.fit = T, interval = "confidence")
fitted$body <- fitted1$fit[,1]
fitted$upper_body <- fitted1$fit[,3]
fitted$lower_body <- fitted1$fit[,2]
head(fitted)


# Mono #
qm_fig1 <- ggplot() +
  geom_boxplot(data = qm, aes(x = clump.size, group = clump.size, y=body), outlier.shape = NA, 
               width = 4, lwd = 0.5) +
  geom_line(data = fitted, aes(x = clump.size, y = body), lwd = 1) +
  geom_ribbon(data = fitted, aes(ymax = upper_body, ymin = lower_body, x = clump.size), alpha = 0.3) +
  xlim(-1,315)+ 
  ylim(0.03,0.2) +
  theme_classic() +
  ggtitle("(a) QM") + xlab("Aggregate size") + ylab("Tissue mass/shell mass (g)") + 
  theme(plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank())
qm_fig1

# Colour #
qm_fig2 <- ggplot() +
  geom_boxplot(data = qm, aes(x = clump.size, group = clump.size, y=body), outlier.shape = NA, 
               width = 4, lwd = 0.5) +
  geom_line(data = fitted, aes(x = clump.size, y = body), col = "coral2", lwd = 1) +
  geom_ribbon(data = fitted, aes(ymax = upper_body, ymin = lower_body, x = clump.size), alpha = 0.3, fill = "coral2") +
  xlim(-1,315)+ 
  ylim(0.03,0.2) +
  theme_classic() +
  ggtitle("(a) QM") + xlab("Aggregate size") + ylab("Tissue mass/shell mass (g)") + 
  theme(plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank())
qm_fig2

sum(qm$clump.size == 1)
sum(qm$clump.size == 28)
sum(qm$clump.size == 38)
sum(qm$clump.size == 48)
sum(qm$clump.size == 92)
sum(qm$clump.size == 99)
sum(qm$clump.size == 116)
sum(qm$clump.size == 282)
sum(qm$clump.size == 310)

######################
#### Zebra mussel ####
######################

zm <- as.data.frame(read.csv("../data/standardized_ZM.csv", header = TRUE, stringsAsFactors = F))
head(zm)
zm$clump.size <- as.numeric(zm$clump.size)
sum(zm$clump.size == 1)

# Fitting segmented model to linear model 
zm_model1 <- lm(body~clump.size*size, zm)
par(mfrow = c(2, 2))
plot(zm_model1)
summary(zm_model1)
anova(zm_model1)

zm_model2 <- segmented(zm_model1, seg.Z = ~clump.size)
summary(zm_model2)
davies.test(zm_model1, seg.Z=~clump.size)

zm_model3 <- lm(body~clump.size+size, zm)
par(mfrow = c(2, 2))
plot(zm_model3)
summary(zm_model3)
anova(zm_model3)

zm_model4 <- segmented(zm_model3, seg.Z = ~clump.size)
summary(zm_model4)
davies.test(zm_model4, seg.Z=~clump.size)

zm_model <- lm(size ~ clump.size, zm)
anova(zm_model)

## Plotting ##
# Predicted values for constant shell length 
fitted3 <- data.frame(clump.size = seq(1, 310, 1))
fitted3$size <- rep(mean(zm$size), nrow(fitted))
head(fitted3)

fitted4 <- predict(zm_model3, newdata = fitted3, se.fit = T, interval = "confidence")
fitted3$body <- fitted4$fit[,1]
fitted3$upper_body <- fitted4$fit[,3]
fitted3$lower_body <- fitted4$fit[,2]
head(fitted3)

# Mono #
zm_fig1 <- ggplot(data=zm, aes(x=clump.size, y=body)) +
  geom_boxplot(data=zm, aes(x=clump.size, y=body, group=clump.size), outlier.shape = NA, width=4, lwd=0.5) +
  guides() +
  geom_line(data = fitted3, aes(x = clump.size, y = body), lwd = 1) +
  geom_ribbon(data = fitted3, aes(ymax = upper_body, ymin = lower_body, x = clump.size), alpha = 0.3) +
  scale_fill_grey(start=0.5, end=0.8) +
  xlim(-1,315)+ ylim(0.03,0.18) +
  theme_classic() +
  ggtitle("(b) ZM") + xlab("Aggregate size") + ylab("Tissue mass/shell mass (g)") + 
  theme(plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank())
zm_fig1


pdf("../results/mono/Clump_size_mono.pdf", width = 7, height = 10)
grid.arrange(qm_fig1,zm_fig1)
dev.off()


# Colour #
zm_fig2 <- ggplot(data=zm, aes(x=clump.size, y=body)) +
  geom_boxplot(data=zm, aes(x=clump.size, y=body, group=clump.size), outlier.shape = NA, width=4, lwd=0.5) +
  guides() +
  geom_line(data = fitted3, aes(x = clump.size, y = body), col = "chartreuse3", lwd = 1) +
  geom_ribbon(data = fitted3, aes(ymax = upper_body, ymin = lower_body, x = clump.size), alpha = 0.3, fill = "chartreuse3") +
  scale_fill_grey(start=0.5, end=0.8) + xlim(-1,315) + ylim(0.03,0.18) +
  theme_classic() +
  ggtitle("(b) ZM") + xlab("Aggregate size") + ylab("Tissue mass/shell mass (g)") + 
  theme(plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank())
zm_fig2

pdf("../results/colour/Clump_size_colour.pdf", width = 7, height = 10)
grid.arrange(qm_fig2, zm_fig2)
dev.off()

sum(zm$clump.size == 1)
sum(zm$clump.size == 28)
sum(zm$clump.size == 38)
sum(zm$clump.size == 48)
sum(zm$clump.size == 92)
sum(zm$clump.size == 99)
sum(zm$clump.size == 116)
sum(zm$clump.size == 282)
sum(zm$clump.size == 310)


### Species ratio vs clump size ###
# Overall ratio
sum(qm$clump.size != 1)/(sum(qm$clump.size != 1) + sum(zm$clump.size != 1))

# Ratio per clump size 
ratio <- as.data.frame(read.csv("../data/Density vs ratio.csv", header = TRUE, stringsAsFactors = F))
ratio

ratio_model <- lm(ratio~no, data = ratio)
summary(ratio_model)
anova(ratio_model)

