###########################
##### Data wrangling ######
###########################
rm(list = ls())

########### QM ############
## Remove unnecessary data from main file "QM dry mass" ##
qm <- as.data.frame(read.csv("../data/raw/QM dry mass.csv", header = T, stringsAsFactors=F))
head(qm)

qm <- qm[,-c(2, 4, 6, 7, 9, 10, 11, 12, 14, 15, 16, 17, 19)]
head(qm)

colnames(qm)[3] <- "tissue"
colnames(qm)[4] <- "shell"
colnames(qm)[5] <- "clump.size"
colnames(qm)[6] <- "type"
head(qm)

qm$clump.size[qm$type == "Singleton"] = 1
qm$clump.size[qm$clump.size == 222] = 310


## Add data of layer position from "QM layer" ##
qm_layer <- as.data.frame(read.csv("../data/raw/QM dry mass layers.csv", header = T, stringsAsFactors = F))
head(qm_layer)

qm$layer <- NA

qm_bottom <- c() # Making a list of ID of individuals in the Bottom layer
for(i in 1:nrow(qm_layer)){
  if(qm_layer$Layer[i] == "Bottom"){
    qm_bottom <- c(qm_bottom, qm_layer$ID[i])
  }
}

qm_top <- c() # Making a list of ID of individuals in the Top layer
for(i in 1:nrow(qm_layer)){
  if(qm_layer$Layer[i] == "Top"){
    qm_top <- c(qm_top, qm_layer$ID[i])
  }
}

# Label individuals in main data as Bottom in layer column from list of ID 
for(i in qm_bottom){     
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$layer[j] <- "Bottom"
    }
  }
}

# Label individuals in main data as Top in layer column from list of ID 
for(i in qm_top){
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$layer[j] <- "Top"
    }
  }
}
head(qm)

## Add data of fouling status from "QM fouling" ##
qm_fouling <- as.data.frame(read.csv("../data/raw/Fouling QM.csv", header = T, stringsAsFactors = F))
head(qm_fouling)
qm$fouling <- NA

# Making a list of ID of individuals with different fouling status 
qm_QM.f <- c()
for(i in 1:nrow(qm_fouling)){
  if(qm_fouling$situ[i] == "QM.f"){
    qm_QM.f <- c(qm_QM.f, qm_fouling$ID[i])
  }
}

qm_ZM.f <- c()
for(i in 1:nrow(qm_fouling)){
  if(qm_fouling$situ[i] == "ZM.f"){
    qm_ZM.f <- c(qm_ZM.f, qm_fouling$ID[i])
  }
}

qm_f.QM <- c()
for(i in 1:nrow(qm_fouling)){
  if(qm_fouling$situ[i] == "f.QM"){
    qm_f.QM <- c(qm_f.QM, qm_fouling$ID[i])
  }
}

qm_f.ZM <- c()
for(i in 1:nrow(qm_fouling)){
  if(qm_fouling$situ[i] == "f.ZM"){
    qm_f.ZM <- c(qm_f.ZM, qm_fouling$ID[i])
  }
}

# Label individuals in the main data fouling column following the list of ID 
for(i in qm_QM.f){
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$fouling[j] <- "Fouled_by_QM"
    }
  }
}

for(i in qm_ZM.f){
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$fouling[j] <- "Fouled_by_ZM"
    }
  }
}

for(i in qm_f.QM){
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$fouling[j] <- "Fouling_QM"
    }
  }
}


for(i in qm_f.ZM){
  for(j in 1:nrow(qm)){
    if(qm$ID[j] == i){
      qm$fouling[j] <- "Fouling_ZM"
    }
  }
}

head(qm)


## Add data from preliminary study ##
prelim <- as.data.frame(read.csv("../data/raw/Sports site.csv", header = T, stringsAsFactors = F))
head(prelim)

prelim <- prelim[, -c(2,3,5,8:22,24,26:28)] 
head(prelim)

colnames(prelim)[4] <- "size" 
colnames(prelim)[5] <- "tissue"
colnames(prelim)[6] <- "shell"
prelim <- subset(prelim, Site == "s1r2")
prelim$clump.size <- nrow(prelim)
head(prelim)

pre_qm <- subset(prelim, Species == "QM")
nrow(pre_qm)
pre_qm <- pre_qm[, -c(2,3)]
pre_qm$type <- rep("Mono", nrow(pre_qm))
pre_qm$layer <- rep(NA, nrow(pre_qm))
pre_qm$fouling <- rep(NA, nrow(pre_qm))
head(pre_qm)

qm <- rbind(qm, pre_qm)
head(qm)
nrow(qm)

write.csv(qm, "../data/raw_QM.csv", row.names = F) 



#### ZM ####
# Load main file "ZM dry mass"
data <- as.data.frame(read.csv("../data/raw/ZM dry mass Sport.csv", header = T, stringsAsFactors=F))
head(data)

data <- data[,-c(2, 7, 8, 10, 11, 12, 13, 14)]
head(data)

colnames(data)[6] <- "clump.size"
colnames(data)[5] <- "type"
head(data)

data$type[data$type == "sing"] = "Singleton"
data$type[data$type == "mono"] = "Mono"
data$type[data$type == "clump"] = "Clump"
data$clump.size[data$type == "Singleton"] = 1
data$clump.size[data$clump.size == 222] = 310


data <- data[,c("ID", "size", "tissue", "shell", "clump.size", "type")]
head(data)

# Add data of layer position from "ZM layer"
layer <- as.data.frame(read.csv("../data/raw/ZM dry mass layer.csv", header = T, stringsAsFactors = F))
head(layer)

data$layer <- NA
head(data)

bottom <- c()
for(i in 1:nrow(layer)){
  if(layer$Layer[i] == "Bottom"){
    bottom <- c(bottom, layer$ID[i])
  }
}

top <- c()
for(i in 1:nrow(layer)){
  if(layer$Layer[i] == "Top"){
    top <- c(top, layer$ID[i])
  }
}

for(i in bottom){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$layer[j] <- "Bottom"
    }
  }
}

for(i in top){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$layer[j] <- "Top"
    }
  }
}
head(data)

# Add data of fouling status from "ZM fouling"
fouling <- as.data.frame(read.csv("../data/raw/ZM fouling.csv", header = T, stringsAsFactors = F))
head(fouling)

data$fouling <- NA

QM.f <- c()
for(i in 1:nrow(fouling)){
  if(fouling$situ[i] == "QM.f"){
    QM.f <- c(QM.f, fouling$ID[i])
  }
}

ZM.f <- c()
for(i in 1:nrow(fouling)){
  if(fouling$situ[i] == "ZM.f"){
    ZM.f <- c(ZM.f, fouling$ID[i])
  }
}

f.QM <- c()
for(i in 1:nrow(fouling)){
  if(fouling$situ[i] == "f.QM"){
    f.QM <- c(f.QM, fouling$ID[i])
  }
}

f.ZM <- c()
for(i in 1:nrow(fouling)){
  if(fouling$situ[i] == "f.ZM"){
    f.ZM <- c(f.ZM, fouling$ID[i])
  }
}

for(i in QM.f){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$fouling[j] <- "Fouled_by_QM"
    }
  }
}

for(i in ZM.f){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$fouling[j] <- "Fouled_by_ZM"
    }
  }
}

for(i in f.QM){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$fouling[j] <- "Fouling_QM"
    }
  }
}


for(i in f.ZM){
  for(j in 1:nrow(data)){
    if(data$ID[j] == i){
      data$fouling[j] <- "Fouling_ZM"
    }
  }
}
head(data)

## Add data from preliminary study ##
pre_zm <- subset(prelim, Species == "ZM")
nrow(pre_zm)
pre_zm <- pre_zm[, -c(2,3)]
head(pre_zm)
pre_zm$type <- rep("Mono", nrow(pre_zm))
pre_zm$layer <- rep(NA, nrow(pre_zm))
pre_zm$fouling <- rep(NA, nrow(pre_zm))
head(pre_zm)

data <- rbind(data, pre_zm)
head(data)
nrow(data)

write.csv(data, "../data/raw_ZM.csv", row.names = F) 

