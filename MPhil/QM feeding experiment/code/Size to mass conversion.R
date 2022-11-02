data <- read.csv(file.choose(), header = TRUE)
data
eq1 <- lm(tissue~size,data=data)
summary(eq1)
plot(tissue~size,data=data)
library(ggplot2)
ggplot(data, aes(x=size, y=tissue)) + 
  geom_smooth(method=lm, se=TRUE, colour="black", size=1) 

eq2 <- lm(tissue~wetmass, data=data)
summary(eq2)
plot(tissue~wetmass, data=data)
ggplot(data, aes(x=wetmass, y=tissue)) + 
  geom_smooth(method=lm, se=TRUE, colour="black", size=1)

std <- function(x) sd(x)/sqrt(length(x))
std(data$size)
mean(data$size)
std(data$tissue)
mean(data$tissue)
