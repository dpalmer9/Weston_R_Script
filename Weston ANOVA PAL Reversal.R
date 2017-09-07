## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
##library(nlm3)

## Acquire data ##
##raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\PD Reversal.csv')
raw.data = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\dPAL.csv')
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
    corr.lat.list = dataset[a,c(204:239)]
    colnames(corr.lat.list) = NULL
    corr.lat.list = as.vector(corr.lat.list[1, ])
    corr.lat.list = corr.lat.list[!is.na(corr.lat.list)]
    rew.lat.list = dataset[a,c(280:315)]
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
raw.dataclean = raw.data[ ,c(2:7,11,14:17,318,319)]



## Separate By Strain ##
data.3x.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 3XTG") & (raw.dataclean$Age..months.=="4")), ]
data.5x.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 5FAD") & (raw.dataclean$Age..months.=="4")), ]
data.app.4m = raw.dataclean[which((raw.dataclean$Mouse.strain==" APPPS1") & (raw.dataclean$Age..months.=="4")), ]

data.3x.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 3XTG") & (raw.dataclean$Age..months.=="10")), ]
data.5x.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" 5FAD") & (raw.dataclean$Age..months.=="10")), ]
data.app.10m = raw.dataclean[which((raw.dataclean$Mouse.strain==" APPPS1") & (raw.dataclean$Age..months.=="10")), ]
## Create Object Sets for Separate Measures ##
Data.GenerateMeasureList.Function <- function(dataset){
  datalist = list()
  
  datalist$totaltime = dataset[ ,c(1,4,5,7,8)]
  datalist$totaltrial = dataset[ ,c(1,4,5,7,9)]
  datalist$corrections = dataset[ ,c(1,4,5,7,10)]
  datalist$acc = dataset[ ,c(1,4,5,7,11)]
  datalist$corrlat = dataset[ ,c(1,4,5,7,12)]
  datalist$rewlat = dataset[ ,c(1,4,5,7,13)]
  
  return(datalist)
}

data.3x.4m.list = Data.GenerateMeasureList.Function(data.3x.4m)
data.5x.4m.list = Data.GenerateMeasureList.Function(data.5x.4m)
data.app.4m.list = Data.GenerateMeasureList.Function(data.app.4m)

data.3x.10m.list = Data.GenerateMeasureList.Function(data.3x.10m)
data.5x.10m.list = Data.GenerateMeasureList.Function(data.5x.10m)
data.app.10m.list = Data.GenerateMeasureList.Function(data.app.10m)
## Transform Data Long to Wide - Listwise Mode ##
Data.ReshapeLongtoWide.Function = function(dataset){
  colnames(dataset) = c('AnimalID','Genotype','Gender','Week','Measure')
  dataset = data.frame(dataset, stringsAsFactors = FALSE)
  data.melt = melt(dataset, id.vars = c('AnimalID','Genotype','Gender','Week'))
  data.cast = dcast(data.melt, AnimalID + Genotype + Gender ~ variable + Week, fun.aggregate = mean, na.rm=TRUE)
}

data.3x.4m.list = lapply(data.3x.4m.list, Data.ReshapeLongtoWide.Function)
data.5x.4m.list = lapply(data.5x.4m.list, Data.ReshapeLongtoWide.Function)
data.app.4m.list = lapply(data.app.4m.list, Data.ReshapeLongtoWide.Function)

data.3x.10m.list = lapply(data.3x.10m.list, Data.ReshapeLongtoWide.Function)
data.5x.10m.list = lapply(data.5x.10m.list, Data.ReshapeLongtoWide.Function)
data.app.10m.list = lapply(data.app.10m.list, Data.ReshapeLongtoWide.Function)

## Generate iData ##
idata = as.data.frame(matrix(c(1:9)))
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
data.3x.10m.anova = lapply(data.3x.10m.list, Data.GenerateLMANOVA.Function, idata=idata)

data.5x.4m.anova = lapply(data.5x.4m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.5x.10m.anova = lapply(data.5x.10m.list, Data.GenerateLMANOVA.Function, idata=idata)

data.app.4m.anova = lapply(data.app.4m.list, Data.GenerateLMANOVA.Function, idata=idata)
data.app.10m.anova = lapply(data.app.10m.list, Data.GenerateLMANOVA.Function, idata=idata)


## SME ##

SME.Function.PAL = function(dataset, idata, measures, within.levels){
  SME.list = list()
  repeated.measure = TRUE
  dataset = dataset[complete.cases(dataset), ]
  idesign = "~Day"
  if((is.element('Day',measures) == FALSE) & (is.element('Day',within.levels) == FALSE)){
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=4))
    new.data.col[ ,c(1:3)] = dataset[ ,c(1:3)]
    new.data.col[ ,4] = apply(dataset[ ,c(4:12)], 1, mean, na.rm = TRUE)
    dataset = new.data.col
    colnames(dataset) = c('AnimalID','Genotype','Gender','Measure')
    idesign = NA
    repeated.measure = FALSE
  }
  if((is.element('Day',measures) == TRUE) & (is.element('Day',within.levels) == FALSE)){
    idesign = "~Day"
    idesign = as.formula(idesign)
    idata.edit = idata
    repeated.measure = TRUE
  }
  if((is.element('Day',measures) == FALSE) & (is.element('Day',within.levels) == TRUE)){
    idesign = NA
    idata.edit = idata
    repeated.measure = FALSE
  }
  if(isTRUE(within.levels == "Gender")){
    filter.list = list(' M'," F")
    filter.col = 3
  }else if(isTRUE(within.levels == "Genotype")){
    filter.list = list(' w',' t')
    filter.col = 2
  }else if(isTRUE(within.levels == "Day")){
    filter.list = list(4,5,6,7,8,9,10,11,12)
    idata.edit = NA
    idesign= NA
    filter.col = NA
  }
  if(is.element('Day',measures) == TRUE){
    measures = measures[measures != "Day"]
  }
  measure.string = "data.matrix ~"
  if(length(measures) == 0){
    measure.string = "data.matrix ~ 1"
  }else{
    for(a in 1:length(measures)){
      if((measure.string == "data.matrix ~") == TRUE){
        measure.string = paste(measure.string,measures[a],sep=" ")
      }else{
        measure.string = paste(measure.string,measures[a],sep=" * ")
      }
    }
  }
  measure.string = as.formula(measure.string)
  if((is.na(filter.col) == TRUE) & (repeated.measure == TRUE)){
    for(a in 1:length(filter.list)){
      current.filter = unlist(filter.list[a])
      dataset.anova = dataset[ ,c(1:3, current.filter)]
      dataset.anova = dataset.anova[complete.cases(dataset.anova), ]
      data.matrix = as.matrix(dataset.anova[ ,4:ncol(dataset.anova)])
      data.lm = lm(measure.string, data = dataset.anova)
      data.anova = Anova(data.lm, idata=idata.edit,idesign= ~Day, type="III")
      data.summary = summary(data.anova, multivariate=FALSE)
      #data.summary = data.summary$univariate.tests
      SME.list[[a]] = data.summary
    }
  }else if((is.na(filter.col) == FALSE) & (repeated.measure == TRUE)){
    for(a in 1:length(filter.list)){
      dataset.anova = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ]
      dataset.anova = dataset.anova[complete.cases(dataset.anova), ]
      data.matrix = as.matrix(dataset.anova[, 4:ncol(dataset.anova)])
      data.lm = lm(measure.string, data = dataset.anova)
      data.anova = Anova(data.lm, idata=idata.edit,idesign=~Day, type="III")
      data.summary = summary(data.anova, multivariate=FALSE)
      #data.summary = data.summary$univariate.tests
      SME.list[[a]] = data.summary
    }
  }else if((is.na(filter.col) == TRUE) & (repeated.measure == FALSE)){
    for(a in 1:length(filter.list)){
      data.matrix = as.matrix(dataset[ ,unlist(filter.list[a])])
      data.lm = lm(measure.string, data = dataset)
      data.anova = Anova(data.lm, type="III")
      #data.summary = summary(data.anova, multivariate=FALSE)
      SME.list[[a]] = data.anova
    }
  }else if((is.na(filter.col) == FALSE) & (repeated.measure == FALSE)){
    for(a in 1:length(filter.list)){
      data.matrix = as.matrix(dataset[which(dataset[filter.col] == unlist(filter.list[a])), 4])
      data.lm = lm(measure.string, data = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ])
      data.anova = Anova(data.lm, type="III")
      SME.list[[a]] = data.anova
    }
  }
  return(SME.list)
}

SME.Function.PAL(data.3x.10m.list$corrections , idata = idata, measures= c('Gender'), within.levels = "Day")


SME.Function.PAL(test, idata = idata, measures= c('Genotype'), within.levels = "Day")
