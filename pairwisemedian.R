load("~/Documents/project/subtypedata.Rdata")

medianHER2Normaldata <- numeric(length = length(HER2Normaldata[,1]))
medianlumAlumBdata <- numeric(length = length(lumAlumBdata[,1]))
medianlumAHER2data <- numeric(length = length(lumAHER2data[,1]))
medianHER2Basaldata <- numeric(length = length(HER2Basaldata[,1]))
medianlumBNormaldata <- numeric(length = length(lumBNormaldata[,1]))
medianlumANormaldata <- numeric(length = length(lumANormaldata[,1]))
medianlumABasaldata <- numeric(length = length(lumABasaldata[,1]))
medianlumBHER2data <- numeric(length = length(lumBHER2data[,1]))
medianBasalNormaldata <- numeric(length = length(BasalNormaldata[,1]))
medianlumBBasaldata <- numeric(length = length(lumBBasaldata[,1]))

for (i in 1: nrow(lumABasaldata)){ 
  medianlumABasaldata[i] <-median(as.matrix(lumABasaldata[i,]))
  medianBasalNormaldata[i] <-median(as.matrix(BasalNormaldata[i,]))
  medianlumAHER2data[i] <-median(as.matrix(lumAHER2data[i,]))
  medianlumAlumBdata[i] <-median(as.matrix(lumAlumBdata[i,]))
  medianlumANormaldata[i] <-median(as.matrix(lumANormaldata[i,]))
  medianlumBNormaldata[i] <-median(as.matrix(lumBNormaldata[i,]))
  medianlumBHER2data[i] <-median(as.matrix(lumBHER2data[i,]))
  medianlumBBasaldata[i] <-median(as.matrix(lumBBasaldata[i,]))
  medianHER2Basaldata[i] <-median(as.matrix(HER2Basaldata[i,]))
  medianHER2Normaldata[i] <-median(as.matrix(HER2Normaldata[i,]))
  }

save(medianlumAlumBdata,medianlumAHER2data, medianlumABasaldata, medianlumANormaldata, medianlumBNormaldata, medianlumBHER2data,medianlumBBasaldata, 
medianHER2Basaldata, medianHER2Normaldata, medianBasalNormaldata, file="~/Documents/project/pairwisemediandata.Rdata")
