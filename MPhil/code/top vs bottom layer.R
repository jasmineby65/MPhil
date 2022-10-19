###################################################
####### Body condition vs vertical position #######
###################################################

rm(list = ls())

require(ggplot2)
require(gridExtra)
require(ggpubr)
require(interactions)
require(lsmeans)

#######################
#### Quagga mussel ####
#######################
qm <- as.data.frame(read.csv("../data/standardized_QM.csv", header = TRUE, stringsAsFactors = F))
head(qm)

qm <- subset(qm, !(is.na(qm$layer)))
qm$clump.size <- as.factor(qm$clump.size)

qm_model1 <- lm(body~layer*size*clump.size, data=qm)
summary(qm_model1)
anova(qm_model1)

lstrends(qm_model1, c("layer", "clump.size"), var ="size")

TukeyHSD((aov(qm_model1)), "layer:clump.size")

qm_model2 <- update(qm_model1, ~. -layer:size:clump.size)
summary(qm_model2)
anova(qm_model2, qm_model1)

qm_model3 <- update(qm_model2, ~. -size:clump.size)
summary(qm_model3)
anova(qm_model3)
TukeyHSD((aov(qm_model3)), "layer:clump.size")

qm_model2 <- lm(body~layer*clump.size, data = qm)
cat_plot(qm_model2, pred = layer, modx = clump.size)
summary(qm_model2)
TukeyHSD((aov(qm_model3)), "layer:clump.size")


# Effect of layer on shell length 
qm_model3 <- lm(size ~ layer*clump.size, data = qm)
summary(qm_model3)
anova(qm_model3)

qm_model4 <- lm(size ~ layer+clump.size, data = qm)
summary(qm_model4)

qm_model5 <- lm(size ~ layer, data = qm)
summary(qm_model5)
anova(qm_model5)

# Plotting #
# Interaction effect - colour #
pdf("../results/colour/layer_interaction_colour.pdf", width = 7.5, height = 5)
interact_plot(qm_model1, pred = size, modx = layer, mod2 = clump.size,
              mod2.labels = c("Clump size = 282", "Clump size = 310"), 
              x.label = "Shell length (mm)", y.label = "Tissue mass/shell mass (g)", legend.main = "Position") +
  theme(legend.position = c(0.9, 0.18), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) 
dev.off()


# Interaction effect - mono #
pdf("../results/mono/layer_interaction_mono.pdf", width = 7.5, height = 5)
interact_plot(qm_model1, pred = size, modx = layer, mod2 = clump.size, colors = "Greys",
              mod2.labels = c("Clump size = 282", "Clump size = 320"), x.label = "Shell length (mm)", 
              y.label = "Tissue mass/shell mass (g)", legend.main = "Position") +
  theme(legend.position = c(0.9, 0.18), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) 
dev.off()


# Layer & clump.size vs body condition # 
plot <- data.frame(size = c("282", "282", "310", "310"))
plot$layer <- c("bottom", "top", "bottom", "top")

a <- qm$body[qm$layer == "Bottom" & qm$clump.size == "282"]
b <- qm$body[qm$layer == "Top" & qm$clump.size == "282"]
c <- qm$body[qm$layer == "Bottom" & qm$clump.size == "310"]
d <- qm$body[qm$layer == "Top" & qm$clump.size == "310"]

length(a)
length(b)
length(c)
length(d)

qm_situ <- list(a,b,c,d)

plot$mean <- NA
plot$CI_up <- NA
plot$CI_low <- NA

for(i in 1:4){
  data <- qm_situ[[i]]
  plot$mean[i] <- mean(data)
  
  model <- lm(data~1)
  conf <- confint(model, level = 0.95)
  plot$CI_low[i] <- conf[1]
  plot$CI_up[i] <- conf[2]
}
plot


# Colour
qm_fig1 <-ggplot(data=plot, aes(x=size , y=mean, colour = layer)) +
  geom_point(position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  theme_classic() +
  xlab("Aggregate size")+
  ylab("tissue mass/shell mass (g)") + ylim(0.068, 0.102) +
  ggtitle("(a) QM") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(),  
        legend.title = element_blank(), 
        legend.position = c(0.7, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(2, "cm"))
qm_fig1

# Mono
qm_fig2 <-ggplot(data=plot, aes(x=size , y=mean, linetype = layer, colour = layer)) +
  geom_point(position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_colour_grey(start=0, end = 0) +
  theme_classic() +
  xlab("Aggregate size")+
  ylab("tissue mass/shell mass (g)") + ylim(0.068, 0.102) +
  ggtitle("(a) QM") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(),
        legend.title = element_blank(), 
        legend.position = c(0.7, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(2, "cm"))
qm_fig2



######################
#### Zebra mussel ####
######################
zm <- as.data.frame(read.csv("../data/standardized_ZM.csv", header = TRUE, stringsAsFactors = F))
head(zm)

zm <- subset(zm, !(is.na(zm$layer)))
zm$clump.size <- as.factor(zm$clump.size)

zm_model1 <- lm(body~layer*size*clump.size, data=zm)
summary(zm_model1)
anova(zm_model1)

zm_model2 <- lm(body~layer+clump.size+size, data=zm)
summary(zm_model2)
anova(zm_model2)

zm_model3 <- lm(body~layer, data=zm)
summary(zm_model3)
anova(zm_model3)


# Effect of layer on shell length 
zm_model4 <- lm(size ~ layer*clump.size, data=zm)
summary(zm_model4)
anova(zm_model4)

zm_model5 <- lm(size ~ layer+clump.size, data=zm)
summary(zm_model5)
anova(zm_model5)

zm_model6 <- lm(size ~ layer, data=zm)
anova(zm_model6)

mean(zm$size[zm$layer == "Top"])
sqrt(var(zm$size[zm$layer == "Top"]) / length(zm$size[zm$layer == "Top"]))
mean(zm$size[zm$layer == "Bottom"])
sqrt(var(zm$size[zm$layer == "Bottom"]) / length(zm$size[zm$layer == "Bottom"]))

# Plotting 

zm_plot <- data.frame(size = c("282", "282", "310", "310"))
zm_plot$layer <- c("bottom", "top", "bottom", "top")

e <- zm$body[zm$layer == "Bottom" & zm$clump.size == "282"]
f <- zm$body[zm$layer == "Top" & zm$clump.size == "282"]
g <- zm$body[zm$layer == "Bottom" & zm$clump.size == "310"]
h <- zm$body[zm$layer == "Top" & zm$clump.size == "310"]
situ <- list(e,f,g,h)

length(e)
length(f)
length(g)
length(h)

zm_plot$mean <- NA
zm_plot$CI_low <- NA
zm_plot$CI_up <- NA

for(i in 1:4){
  data <- situ[[i]]
  zm_plot$mean[i] <- mean(data)

  model <- lm(data~1)
  conf <- confint(model, level = 0.95)
  zm_plot$CI_low[i] <- conf[1]
  zm_plot$CI_up[i] <- conf[2]
}
zm_plot


# Colour
zm_fig1 <-ggplot(data=zm_plot, aes(x=size , y=mean, colour = layer)) +
  geom_point(position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  theme_classic() +
  xlab("Aggregate size")+
  ylab("tissue mass/shell mass (g)") + ylim(0.045, 0.12) +
  ggtitle("(b) ZM") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(),
        legend.title = element_blank(), 
        legend.position = c(0.7, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(2, "cm"))
zm_fig1

pdf("../results/colour/layer_colour.pdf", width = 6, height = 11)
grid.arrange(qm_fig1, zm_fig1)
dev.off()

# Mono
zm_fig2 <-ggplot(data=zm_plot, aes(x=size , y=mean, linetype = layer, colour = layer)) +
  geom_point(position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_colour_grey(start=0, end = 0) +
  theme_classic() +
  xlab("Aggregate size")+
  ylab("tissue mass/shell mass (g)") + ylim(0.045, 0.12) +
  ggtitle("(b) ZM") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(0.7, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(2, "cm"))
zm_fig2

pdf("../results/mono/layer_mono.pdf", width = 6, height = 11)
grid.arrange(qm_fig2, zm_fig2)
dev.off()
