rm(list = ls())
require(smatr)

##############################
####### Quagga Mussel ########
##############################

QM <- as.data.frame(read.csv("../data/raw_QM.csv", header = T, stringsAsFactors=F))
head(QM)

### Standardizing dry mass ###
## Tissue ##
# ln-transformation
QM$lntissue <- log(QM$tissue) 
QM$lnsize <- log(QM$size) 

# Fitting SMA regression 
qm_sma1 <- sma(lntissue ~ lnsize, data = QM) 
qm_sma1
qm_bsma1 <- qm_sma1$coef[[1]]$`coef(SMA)`[2] # obtaining b_sma value

# Calculating standardized dry mass
qm_size0 <- mean(QM$size) # mean shell length
QM$tissue <- (QM$tissue*((qm_size0/QM$size)^qm_bsma1))
head(QM)


## Shell ## 
# ln-transformation
QM$lnshell <- log(QM$shell) 

# Fitting SMA regression 
qm_sma2 <- sma(lnshell ~ lnsize, data = QM) 
qm_sma2
qm_bsma2 <- qm_sma2$coef[[1]]$`coef(SMA)`[2] # obtaining b_sma value

# Calculating standardized dry mass
QM$shell <- (QM$shell*((qm_size0/QM$size)^qm_bsma2))
head(QM)


## Calculating body condition ##
QM$body <- QM$tissue/QM$shell


## Removing outliers ##
QM <- subset(QM, size != max(QM$size))
QM <- subset(QM, ID != "E8")


## Modelling standardised allometry ##
qm_tissue <- lm(log(QM$tissue) ~ log(QM$size))
qm_shell <- lm(log(QM$shell) ~ log(QM$size))

summary(qm_tissue)
confint(qm_tissue)
summary(qm_shell)
confint(qm_shell)


## Save data ##
head(QM)
QM <- QM[, -c(9:11)]
head(QM)
write.csv(QM, "../data/standardized_QM.csv", row.names = F)



###########################
###### Zebra Mussel #######
###########################

ZM <- as.data.frame(read.csv("../data/raw_ZM.csv", header = T, stringsAsFactors=F))
head(ZM)

### Standardizing dry mass ###
## Tissue ##
# ln-transformation
ZM$lntissue <- log(ZM$tissue) 
ZM$lnsize <- log(ZM$size) 
for(i in 1:nrow(ZM)){ # Remove -infinite values 
  if(ZM$lntissue[i] == -Inf){
    ZM$lntissue[i] <- NA
  }
}

# Fitting SMA regression 
zm_sma1 <- sma(lntissue ~ lnsize, data = ZM) 
zm_sma1
zm_bsma1 <- zm_sma1$coef[[1]]$`coef(SMA)`[2] # obtaining b_sma value

# Calculating standardized dry mass
zm_size0 <- mean(ZM$size) # mean shell length
ZM$tissue <- (ZM$tissue*((zm_size0/ZM$size)^zm_bsma1))
head(ZM)


## Shell ## 
# ln-transformation
ZM$lnshell <- log(ZM$shell) 

# Fitting SMA regression 
zm_sma2 <- sma(lnshell ~ lnsize, data = ZM) 
zm_sma2
zm_bsma2 <- zm_sma2$coef[[1]]$`coef(SMA)`[2] # obtaining b_sma value

# Calculating standardized dry mass
ZM$shell <- (ZM$shell*((zm_size0/ZM$size)^zm_bsma2))
head(ZM)

## Calculating body condition ##
ZM$body <- ZM$tissue/ZM$shell


## Modelling standardised allometry ##
lntissue <- data.frame(tissue = log(ZM$tissue), size = log(ZM$size)) 
lntissue <- subset(lntissue, tissue != -Inf)

zm_tissue <- lm(tissue ~ size, data = lntissue)
zm_shell <- lm(log(ZM$shell) ~ log(ZM$size))

summary(zm_tissue)
confint(zm_tissue)
summary(zm_shell)
confint(zm_shell)


## Save data ##
head(ZM)
ZM <- ZM[, -c(9:11)]
head(ZM)
write.csv(ZM, "../data/standardized_ZM.csv", row.names = F)


##########################
######## Plotting ########
##########################

### SMA regression - colour ###
dev.off()
pdf("../results/colour/sma_plot_colour.pdf", width = 5, height = 9)
layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 1),  
       widths = c(1, 1), respect = T)    

## QM ###
plot(qm_sma1, xlim = c(1.8, 3.1), ylim = c(-7, 0), 
     xlab = "ln(shell length)", ylab = "ln(dry mass)",
     col = "coral2", pch = 1, from = 1.7,  to = 3.2)
plot(qm_sma2, col="chartreuse3", pch = 1, add = T, from = 1.7,  to = 3.2)
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       fill=c("coral2", "chartreuse3"), inset = c(0, 1), bty="n")
title(main = "(a) QM", font.main = 1, adj = 0, line = 0.8)

## ZM ## 
plot(zm_sma1, xlim = c(1.7, 2.9), ylim = c(-7, -1),
     xlab = "ln(shell length)", ylab = "ln(dry mass)",
     col = "coral2", pch = 1, from = 1.6, to = 3)
plot(zm_sma2, col="chartreuse3", pch = 1, add = T, from = 1.6, to = 3)
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       fill=c("coral2", "chartreuse3"), inset = c(0, 1), bty="n")
title(main = "(b) ZM", font.main = 1, adj = 0, line = 0.8)

dev.off()

### Allometric scaling - colour ###
pdf("../results/colour/allometric_plot_colour.pdf", width = 5, height = 9)

layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 1),  
       widths = c(1, 1), respect = T)    


## QM ##
plot(log(tissue) ~ log(size), data = QM, col="coral2", 
     xlab = "ln(shell length)", ylab = "ln(standardized dry mass)",
     ylim = c(-6, -1))
abline(lm(log(QM$tissue)~log(QM$size)), col="coral2")
points(log(shell)~log(size), data = QM, col = "chartreuse3")
abline(lm(log(QM$shell)~log(QM$size)), col="chartreuse3")
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       fill=c("coral2", "chartreuse3"), inset = c(0, 1), bty="n")
title(main = "(a) QM", font.main = 1, adj = 0, line = 0.8)

## ZM ##
lntissue <- log(ZM$tissue)
for(i in lntissue){ # Remove -infinite values 
  if(i == -Inf){
    lntissue[lntissue == i] <- NA
  }
}

plot(log(tissue) ~ log(size), data = ZM, col="coral2", 
     ylim = c(-6.5, -2), xlim = c(1.7, 2.9),
     xlab = "ln(shell length)", ylab = "ln(standardized dry mass)")
abline(lm(lntissue~log(ZM$size)), col="coral2")
points(log(shell)~log(size), data = ZM, col = "chartreuse3")
abline(lm(log(ZM$shell)~log(ZM$size)), col="chartreuse3")
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       fill=c("coral2", "chartreuse3"), inset = c(0, 1), bty="n")
title(main = "(b) ZM", font.main = 1, adj = 0, line = 0.8)

dev.off()


### SMA regression - mono ###
pdf("../results/mono/sma_plot_mono.pdf", width = 5, height = 9)
layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 1),  
       widths = c(1, 1), respect = T)    

## QM ###
plot(qm_sma1, xlim = c(1.8, 3.1), ylim = c(-7, 0),
     xlab = "ln(shell length)", ylab = "ln(dry mass)",
     pch = 1, col = "grey12", lty = 1, from = 1.7,  to = 3.2)
plot(qm_sma2, pch = 2, col = "grey12", lty = 1, add = T, from = 1.7,  to = 3.2)
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       lty=c(1, 1), pch = c(1,2), inset = c(0, 1), bty="n")
title(main = "(a) QM", font.main = 1, adj = 0, line = 0.8)

## ZM ## 
plot(zm_sma1, xlim = c(1.7, 2.9), ylim = c(-7, -1),
     xlab = "ln(shell length)", ylab = "ln(dry mass)",
     pch = 1, lty = 1, col = "grey12", from = 1.6, to = 3)
plot(zm_sma2, pch = 2, lty = 1, add = T, col = "grey12", from = 1.6, to = 3)
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       lty=c(1, 1), pch = c(1,2), inset = c(0, 1), bty="n")
title(main = "(b) ZM", font.main = 1, adj = 0, line = 0.8)

dev.off()

### Allometric scaling - mono ###
pdf("../results/mono/allometric_plot_mono.pdf", width = 5, height = 9)

layout(mat = matrix(c(1, 2), 
                    nrow = 2, 
                    ncol = 1),
       heights = c(1, 1),  
       widths = c(1, 1), respect = T)    

## QM ##
plot(log(tissue) ~ log(size), data = QM, col="grey12", 
     xlim = c(1.8, 3.1), ylim = c(-6, -1), pch = 1,
     xlab = "ln(shell length)", ylab = "ln(standardized dry mass)")
abline(lm(log(QM$tissue)~log(QM$size)), col="grey12")
points(log(shell)~log(size), data = QM, col = "grey12", pch = 2)
abline(lm(log(QM$shell)~log(QM$size)), col="grey12")
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       pch=c(1,2), lty = c(1,1), inset = c(0, 1), bty="n")
title(main = "(a) QM", font.main = 1, adj = 0, line = 0.8)

## ZM ##
plot(log(tissue) ~ log(size), data = ZM, col="grey12", pch = 1, 
     xlim = c(1.7, 2.9), ylim = c(-6.5, -2), 
     xlab = "ln(shell length)", ylab = "ln(standardized dry mass)")
abline(lm(lntissue~log(ZM$size)), col="grey12")
points(log(shell)~log(size), data = ZM, col = "grey12", pch = 2)
abline(lm(log(ZM$shell)~log(ZM$size)), col="grey12")
legend('bottomright',c('Tissue','Shell'), horiz = T, xpd = T, 
       pch=c(1,2), lty = c(1,1), inset = c(0, 1), bty="n")
title(main = "(b) ZM", font.main = 1, adj = 0, line = 0.8)

dev.off()

