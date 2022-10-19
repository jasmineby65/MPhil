################################################
####### Body condition vs aggregate type #######
################################################

rm(list = ls())

require(ggplot2)
require(lsmeans)
require(gridExtra)
require(grid)
require(gridtext)
require(smatr)


#######################
#### Quagga mussel ####
#######################

qm <- as.data.frame(read.csv("../data/standardized_QM.csv", header = TRUE, stringsAsFactors = F))
head(qm)

qm_model1 <- lm(body ~ type*size, data=qm)
summary(qm_model1)
anova(qm_model1)

qm_model2 <- lm(body ~ type+size, data=qm)
anova(qm_model2)
anova(qm_model1, qm_model2)
TukeyHSD((aov(qm_model2)), "type")

qm_model <- lm(size ~ type, data = qm)
summary(qm_model)
anova(qm_model)
TukeyHSD((aov(qm_model)), "type")

# Plotting #
qm$Type <- factor(qm$type,
                  levels = c('Singleton','Mono','Clump'),ordered = TRUE)

# Colour 
# qm_fig1 <-ggplot(data=qm, aes(x=Type,y=body)) +
#   geom_boxplot(data=qm, aes(x=Type, y=body, fill = Type),lwd=0.5) +
#   # scale_fill_grey(start=0.5, end=0.8)+
#   theme_classic() +
#   ggtitle("(a) QM") +  scale_fill_manual(values= wes_palette("GrandBudapest2", n = 3)) +
#   ylim(0,0.31) +
#   theme(legend.position = "none") + 
#   theme(aspect.ratio = 1, 
#         plot.title = element_text(size = 15), axis.title = element_blank(),
#         panel.border = element_rect(size = 0.5, fill = NA),
#         axis.line = element_blank(), plot.margin = unit(rep(0.3, 4), "cm"))
# qm_fig1


# Mono
qm_fig2 <-ggplot(data=qm, aes(x=Type,y=body)) +
  geom_boxplot(data=qm, aes(x=Type, y=body),lwd=0.5) +
  theme_classic() +
  ggtitle("(a) QM") + xlab("Aggregate type") + ylab("Tissue mass/shell mass (g)") +
  ylim(0,0.31) +
  theme(legend.position = "none") + 
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15), 
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank())
qm_fig2

sum(qm$type == "Singleton")
sum(qm$type == "Mono")
sum(qm$type == "Clump")

######################
#### Zebra mussel ####
######################
zm <- as.data.frame(read.csv("../data/standardized_ZM.csv", header = TRUE, stringsAsFactors = F))

zm_model1 <- lm(body ~ type*size, data=zm)
summary(zm_model1)
anova(zm_model1)

zm_model2 <- lm(body ~ type+size, data=zm)
anova(zm_model2)
anova(zm_model1, zm_model2)
TukeyHSD((aov(zm_model2)), "type")

zm_model <- lm(size ~ type, data = zm)
summary(zm_model)
anova(zm_model)

# Plotting #
zm$Type <- factor(zm$type,
                  levels = c('Singleton','Mono','Clump'),ordered = TRUE)

# Colour
# zm_fig1 <-ggplot(data=zm, aes(x=Type,y=body)) +
#   geom_boxplot(data=zm, aes(x=Type,y=body, fill = Type),lwd=0.5) +
#   theme_classic() +
#   ggtitle("(b) ZM") +
#   ylim(0,0.4)+scale_fill_manual(values= wes_palette("GrandBudapest2", n = 3)) +
#   theme(legend.position = "none") + 
#   theme(aspect.ratio = 1, 
#         plot.title = element_text(size = 15), axis.title = element_blank(),
#         panel.border = element_rect(size = 0.5, fill = NA),
#         axis.line = element_blank(), plot.margin = unit(rep(0.3, 4), "cm"))
# zm_fig1
# 
# left <- richtext_grob("Tissue mass/shell mass (g)", rot = 90)
# bottom <- richtext_grob("Aggregate type")
# 
# pdf("../results/colour/Aggregate_type_colour.pdf", width = 5, height = 10)
# grid.arrange(qm_fig1, zm_fig1, left = left, bottom = bottom)
# dev.off()

# Mono
zm_fig2 <-ggplot(data=zm, aes(x=Type,y=body)) +
  geom_boxplot(data=zm, aes(x=Type,y=body),lwd=0.5) +
  theme_classic() +
  ggtitle("(b) ZM") + xlab("Aggregate type") + ylab("Tissue mass/shell mass (g)") +
  ylim(-0.03,0.4) +
  theme(legend.position = "none") + 
  theme(aspect.ratio = 1, 
        plot.title = element_text(size = 15),
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.line = element_blank(), plot.margin = unit(rep(0.3, 4), "cm"))
zm_fig2

sum(zm$type == "Singleton")
sum(zm$type == "Mono")
sum(zm$type == "Clump")

pdf("../results/mono/Aggregate_type_mono.pdf", width = 5, height = 10)
grid.arrange(qm_fig2, zm_fig2)
dev.off()
