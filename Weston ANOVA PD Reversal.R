## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
##library(nlm3)

## Acquire data ##
##raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\PD Reversal.csv')
raw.data = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\PD Reversal.csv')
raw.data$End.Summary...Condition = as.numeric(raw.data$End.Summary...Condition)


## Re Calculate Mean Value for Error ##
Data.LatencyCalc.Function <- function(dataset){
  if((("mean.corrlat" %in% colnames(dataset)) & ("mean.rewlat" %in% colnames(dataset))) == FALSE){
    dataset$mean.corrlat = NA
    dataset$mean.rewlat = NA
    corrlat.colnum = which(colnames(dataset) == "mean.corrlat")
    rewlat.colnum = which(colnames(dataset) == "mean.rewlat")
  }
  for(a in 1:nrow(dataset)){
    corr.lat.list = dataset[a,c(17:46)]
    colnames(corr.lat.list) = NULL
    corr.lat.list = as.vector(corr.lat.list[1, ])
    corr.lat.list = corr.lat.list[!is.na(corr.lat.list)]
    rew.lat.list = dataset[a,c(49:78)]
    colnames(rew.lat.list) = NULL
    rew.lat.list = as.vector(rew.lat.list[1, ])
    rew.lat.list = rew.lat.list[!is.na(rew.lat.list)]
    if(length(corr.lat.list) == 0){
      dataset[a,corrlat.colnum] = NA
    }else if (length(corr.lat.list) == 1){
      dataset[a,corrlat.colnum] = corr.lat.list[1]
    }else if (length(corr.lat.list) > 1){
      dataset[a,corrlat.colnum] = mean(as.vector(corr.lat.list), na.rm = TRUE)
    }
    
    if(length(rew.lat.list) == 0){
      dataset[a,rewlat.colnum] = NA
    }else if (length(rew.lat.list) == 1){
      dataset[a,rewlat.colnum] = rew.lat.list[1]
    }else if (length(rew.lat.list) > 1){
      dataset[a,rewlat.colnum] = mean(as.vector(rew.lat.list), na.rm = TRUE)
    }
  }
  return(dataset)
}

raw.data = Data.LatencyCalc.Function(raw.data)
## Remove Individual Latency ##
raw.dataclean = raw.data[ ,c(2:8,10,13:16,81,82)]



## Separate By Strain ##
data.3x.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 3XTG") & (raw.dataclean$Age..months.=="4")), ]
data.5x.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 5FAD") & (raw.dataclean$Age..months.=="4")), ]
data.app.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" APPPS1") & (raw.dataclean$Age..months.=="4")), ]

data.3x.7m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 3XTG") & (raw.dataclean$Age..months.=="7")), ]
data.5x.7m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 5FAD") & (raw.dataclean$Age..months.=="7")), ]
data.app.7m = raw.dataclean[which((raw.dataclean$Mouse.strain==" APPPS1") & (raw.dataclean$Age..months.=="7")), ]

data.3x.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 3XTG") & (raw.dataclean$Age..months.=="10")), ]
data.5x.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 5FAD") & (raw.dataclean$Age..months.=="10")), ]
data.app.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" APPPS1") & (raw.dataclean$Age..months.=="10")), ]
## Create Object Sets for Separate Measures ##
Data.GenerateMeasureList.Function <- function(dataset){
  datalist = list()
  
  datalist$totaltime = dataset[ ,c(1,4:6,8,9)]
  datalist$totaltrial = dataset[ ,c(1,4:6,8,10)]
  datalist$corrections = dataset[ ,c(1,4:6,8,11)]
  datalist$acc = dataset[ ,c(1,4:6,8,12)]
  datalist$corrlat = dataset[ ,c(1,4:6,8,13)]
  datalist$rewlat = dataset[ ,c(1,4:6,8,14)]
  
  return(datalist)
}

data.3x.4m.list = Data.GenerateMeasureList.Function(data.3x.4m)
data.5x.4m.list = Data.GenerateMeasureList.Function(data.5x.4m)
data.app.4m.list = Data.GenerateMeasureList.Function(data.app.4m)

data.3x.7m.list = Data.GenerateMeasureList.Function(data.3x.7m)
data.5x.7m.list = Data.GenerateMeasureList.Function(data.5x.7m)
data.app.7m.list = Data.GenerateMeasureList.Function(data.app.7m)

data.3x.10m.list = Data.GenerateMeasureList.Function(data.3x.10m)
data.5x.10m.list = Data.GenerateMeasureList.Function(data.5x.10m)
data.app.10m.list = Data.GenerateMeasureList.Function(data.app.10m)
## Transform Data Long to Wide - Listwise Mode ##
Data.ReshapeLongtoWide.Function = function(dataset){
  colnames(dataset) = c('AnimalID','Genotype','Gender','Age','Day','Measure')
  dataset = data.frame(dataset, stringsAsFactors = FALSE)
  dataset$Genotype[dataset$Genotype == " w"] = "w"
  dataset$Genotype[dataset$Genotype == " t"] = "t"
  dataset$Genotype = as.factor(dataset$Genotype)
  #dataset$Gender[dataset$Gender == " F"] = "F"
  dataset$Gender[dataset$Gender == " M"] = "M"
  dataset$Gender = as.factor(dataset$Gender)
  dataset$Age = NULL
  data.melt = melt(dataset, id.vars = c('AnimalID','Genotype','Gender','Day'))
  data.cast = dcast(data.melt, AnimalID + Genotype + Gender ~ variable + Day, fun.aggregate = mean, na.rm=TRUE)
}

data.3x.4m.list = lapply(data.3x.4m.list, Data.ReshapeLongtoWide.Function)
data.5x.4m.list = lapply(data.5x.4m.list, Data.ReshapeLongtoWide.Function)
data.app.4m.list = lapply(data.app.4m.list, Data.ReshapeLongtoWide.Function)

data.3x.7m.list = lapply(data.3x.7m.list, Data.ReshapeLongtoWide.Function)
data.5x.7m.list = lapply(data.5x.7m.list, Data.ReshapeLongtoWide.Function)
data.app.7m.list = lapply(data.app.7m.list, Data.ReshapeLongtoWide.Function)

data.3x.10m.list = lapply(data.3x.10m.list, Data.ReshapeLongtoWide.Function)
data.5x.10m.list = lapply(data.5x.10m.list, Data.ReshapeLongtoWide.Function)
data.app.10m.list = lapply(data.app.10m.list, Data.ReshapeLongtoWide.Function)

## Generate iData ##
idata = as.data.frame(matrix(c(1:10)))
colnames(idata) = "Day"
idata$Day = as.factor(idata$Day)


Data.GenerateLMANOVA.Function <- function(dataset,idata){
  dataset = dataset[complete.cases(dataset), ]
  dataset$AnimalID = NULL
  data.depend = dataset[ ,3:ncol(dataset)]
  data.lm = lm(as.matrix(data.depend) ~ Genotype * Gender, data=dataset, na.action = na.omit)
  data.anova = Anova(data.lm, idata=idata,idesign=~Day, type="III")
  return(data.anova)
}

data.3x.4m.anova = lapply(data.3x.4m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.3x.7m.anova = lapply(data.3x.7m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.3x.10m.anova = lapply(data.3x.10m.list, Data.GenerateLMANOVA.Function, idata=idata)

data.5x.4m.anova = lapply(data.5x.4m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.5x.7m.anova = lapply(data.5x.7m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.5x.10m.anova = lapply(data.5x.10m.list, Data.GenerateLMANOVA.Function, idata=idata)

data.app.4m.anova = lapply(data.app.4m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.app.7m.anova = lapply(data.app.7m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.app.10m.anova = lapply(data.app.10m.list, Data.GenerateLMANOVA.Function, idata=idata)

