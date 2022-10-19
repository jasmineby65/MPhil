#################################
####### Fouling vs fouled #######
#################################

rm(list = ls())

require(ggplot2)
require(lsmeans)
require(interactions)

#######################
#### Quagga mussel ####
#######################
qm <- as.data.frame(read.csv("../data/standardized_QM.csv", header = TRUE, stringsAsFactors = F))
qm <- subset(qm, !(is.na(qm$fouling)))
qm$clump.size <- as.numeric(qm$clump.size)
head(qm)

unique(qm$clump.size)

qm$state <- NA
qm$state[qm$fouling == "Fouling_QM" | qm$fouling == "Fouling_ZM"] = "Fouling"
qm$state[qm$fouling == "Fouled_by_QM" | qm$fouling == "Fouled_by_ZM"] = "Fouled"

qm$species <- NA
qm$species[qm$fouling == "Fouling_QM" | qm$fouling == "Fouled_by_QM"] = "QM"
qm$species[qm$fouling ==  "Fouling_ZM" | qm$fouling == "Fouled_by_ZM"] = "ZM"
head(qm)


# Correlation between fouling state and position in clump 
qm_layer <- subset(qm, !(is.na(qm$layer)))

table(qm_layer$layer, qm_layer$state)
sum(table(qm_layer$layer, qm_layer$state))
chisq.test(qm_layer$layer, qm_layer$state)


# Effect of fouling state + species
qm_model1 <- lm(body ~ state*species*size*clump.size, data=qm)
summary(qm_model1)
anova(qm_model1)

qm_model2 <- lm(body ~ state+species+size+clump.size, data=qm)
summary(qm_model2)
anova(qm_model2)

p.adjust(summary(qm_model2)$coefficients[,4], method="bonferroni")

qm_model <- lm(size ~ state+species+clump.size, data = qm)
anova(qm_model)

# Plotting #
plot <- data.frame(species = c("QM", "QM", "ZM", "ZM"))
plot$state <- c("Fouling", "Fouled", "Fouling", "Fouled")

a <- qm$body[qm$fouling == "Fouling_QM"]
b <- qm$body[qm$fouling == "Fouled_by_QM"]
c <- qm$body[qm$fouling == "Fouling_ZM"]
d <- qm$body[qm$fouling == "Fouled_by_ZM"]

qm_situ <- list(a,b,c,d)

length(a)
length(b)
length(c)
length(d)

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
qm_fig1 <-ggplot(data=plot, aes(x=state , y=mean, colour = species)) +
  geom_point(aes(shape = species), position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  theme_classic()  +
  ylim(0.05, 0.115) +
  ggtitle("(a) QM") + 
  xlab("Attachment status") + ylab("Tissue mass/shell mass (g)") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(), 
        legend.title = element_blank(), 
        legend.position = c(0.7, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(2, "cm"))
qm_fig1

# Mono
qm_fig2 <-ggplot(data=plot, aes(x=state , y=mean, linetype = species, colour = species)) +
  geom_point(aes(shape = species), position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, position=position_dodge(width=0.5)) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_colour_grey(start=0, end = 0) +
  theme_classic()  + 
  ylim(0.05, 0.115) +
  ggtitle("(a) QM") + 
  xlab("Attachment status") + ylab("Tissue mass/shell mass (g)") +
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
zm <- subset(zm, !(is.na(zm$fouling)))
zm$clump.size <- as.numeric(zm$clump.size)
head(zm)

zm$species <- NA
zm$species[zm$fouling == "Fouling_QM" | zm$fouling == "Fouled_by_QM"] = "QM"
zm$species[zm$fouling ==  "Fouling_ZM" | zm$fouling == "Fouled_by_ZM"] = "ZM"

zm$state <- NA
zm$state[zm$fouling == "Fouling_QM" | zm$fouling == "Fouling_ZM"] = "Fouling"
zm$state[zm$fouling == "Fouled_by_QM" | zm$fouling == "Fouled_by_ZM"] = "Fouled"

head(zm)


# Correlation between fouling state and position in clump 
zm_layer <- subset(zm, !(is.na(zm$layer)))

table(zm_layer$layer, zm_layer$interaction)
chisq.test(zm_layer$layer, zm_layer$interaction)


# Effect of attachment status
zm_model3 <- lm(body ~ state * size * clump.size, data = zm)
summary(zm_model3)
anova(zm_model3)

zm_model4 <- update(zm_model3, ~. -state:size:clump.size)
summary(zm_model4)

zm_model5 <- update(zm_model4, ~. -size:clump.size)
summary(zm_model5)

zm_model6 <- update(zm_model5, ~. -state:clump.size)
summary(zm_model6)
anova(zm_model6)

zm_model7 <- update(zm_model6, ~. -state:size)
summary(zm_model7)
anova(zm_model7)

lstrends(zm_model3, "state", var = "size")

pdf("../results/colour/zm_interaction_colour.pdf", width = 5, height = 5)
interact_plot(zm_model3, pred = size, modx = state, plot.points = T, data = subset(zm, body != 0), point.shape = T,
              x.label = "Shell length", y.label = "tissue mass/shell mass (g)", legend.main = "Attachment status") +
  theme(legend.position = c(0.8, 0.85), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) + guides(shape = "none")
dev.off()

pdf("../results/mono/zm_interaction_mono.pdf", width = 5, height = 5)
interact_plot(zm_model3, pred = size, modx = state, plot.points = T, data = subset(zm, body != 0), 
              colors = c("black", "black"), point.shape = c(1,2),
              x.label = "Shell length", y.label = "tissue mass/shell mass (g)", legend.main = "Attachment status") +
  theme(legend.position = c(0.8, 0.85), legend.background = element_blank(),
        panel.border = element_rect(size = 0.5, fill = NA), axis.ticks = element_line(size = 0.5),
        panel.grid.major = element_blank()) + guides(shape = "none")
dev.off()

sum(zm$interaction == "Fouled")
sum(zm$interaction == "Fouling")

zm_model2 <- lm(size ~ state + clump.size, data = zm)
anova(zm_model)


# Effect of attached speices
zm_model3 <- lm(body ~ species + size + clump.size, data = zm)
anova(zm_model3)

zm_model4 <- lm(size ~ species + clump.size, data = zm)
anova(zm_model4)


# Plotting #
zm_plot <- data.frame(species = c("QM", "QM", rep("ZM",4)))
zm_plot$state <- c("Fouling", "Fouled", "Fouling","Fouling", "Fouled", "Fouled")

e <- zm$body[zm$fouling == "Fouling_QM"]
f <- zm$body[zm$fouling == "Fouled_by_QM"]
g <- zm$body[zm$fouling == "Fouling_ZM"]
h <- zm$body[zm$fouling == "Fouled_by_ZM"]

zm_situ <- list(e,f,g,h)
zm_situ

length(e)
length(f)
length(g)
length(h)

zm_plot$mean <- NA
zm_plot$CI_up <- NA
zm_plot$CI_low <- NA

for(i in 1:2){
  data <- zm_situ[[i]]
  zm_plot$mean[i] <- mean(data)
  
  model <- lm(data~1)
  conf <- confint(model, level = 0.95)
  zm_plot$CI_low[i] <- conf[1]
  zm_plot$CI_up[i] <- conf[2]
}

zm_plot$mean[3] <- g[1]
zm_plot$mean[4] <- g[2]
zm_plot$mean[5] <- h[1]
zm_plot$mean[6] <- h[2]

zm_plot


# Colour
zm_fig1 <-ggplot(data=zm_plot, aes(x=state , y=mean, linetype = species, colour = species)) +
  geom_point(aes(shape = species), position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, 
                position=position_dodge(width=0.5), show.legend = F) +
  scale_linetype_manual(values = c(1, 2)) +
  theme_classic()  + 
  ylim(0.065, 0.12) +
  ggtitle("(a) ZM") +
  xlab("Attachment status") + ylab("Tissue mass/shell mass (g)") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(),
        legend.title = element_blank(), 
        legend.position = c(0.8, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(1, "cm"))
zm_fig1


pdf("../results/colour/attachment_colour.pdf", width = 6, height = 11)
grid.arrange(qm_fig1, zm_fig1)
dev.off()


# Mono
zm_fig2 <-ggplot(data=zm_plot, aes(x=state , y=mean, linetype = species, colour = species)) +
  geom_point(aes(shape = species), position=position_dodge(width=0.5), size = 3)  +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_up), width = 0.3, lwd = 1, 
                position=position_dodge(width=0.5), show.legend = F) +
  scale_linetype_manual(values = c(1, 2)) +
  scale_colour_grey(start=0, end = 0) +
  theme_classic()  + 
  ylim(0.065, 0.12) +
  ggtitle("(a) ZM") +
  xlab("Attachment status") + ylab("Tissue mass/shell mass (g)") +
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(), 
        legend.title = element_blank(), 
        legend.position = c(0.8, 1.04), legend.direction = "horizontal",
        legend.background = element_blank(), legend.key.size = unit(1, "cm"))
zm_fig2

pdf("../results/mono/attachment_mono.pdf", width = 6, height = 11)
grid.arrange(qm_fig2, zm_fig2)
dev.off()
