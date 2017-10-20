#################################################
#                                               #
#                                               #
#         Weston Anova Script                   #
#             5-CSRTT                           #
#                                               #
#                                               #
#################################################

## Library ##
library(tidyverse)
library(reshape2)
library(car)

## Functions List ##

# Separate Into Strains for Analysis #
Strain.Separation.Function = function(dataset,long.form=0){
  new.data = list()
  if(long.form == 0){
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), ]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="5xFAD"), ]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), ] 
  }else if(long.form == 1){
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), c(2:9,14:19,70,122)]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="5xFAD"), c(2:9,14:19,70,122)]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), c(2:9,14:19,70,122)]
  }
  return(new.data)
}

# Separate Into Separate Measures for Each File + Transform # #
Measure.Separation.Function = function(dataset, long.form=0){
  new.data = list()
  if(long.form == 0){
    new.data$APP$Sessions = dataset$APP[ ,c(2:9)]
    new.data$TG5x$Sessions = dataset$TG5x[ ,c(2:9)]
    new.data$TG3x$Sessions = dataset$TG3x[ ,c(2:9)]
  }else if(long.form == 1){
    new.data$APP$TotalTime = dataset$APP[ ,c(2:8,9)]
    new.data$APP$TotalTrials = dataset$APP[ ,c(2:8,10)]
    new.data$APP$Accuracy = dataset$APP[ ,c(2:8,11)]
    new.data$APP$Omission = dataset$APP[ ,c(2:8,12)]
    new.data$APP$Premature = dataset$APP[ ,c(2:8,13)]
    new.data$APP$Perseverative = dataset$APP[ ,c(2:8,14)]
    new.data$APP$RewardLat = dataset$APP[ ,c(2:8,15)]
    new.data$APP$CorrectLat = dataset$APP[ ,c(2:8,16)]
    
    new.data$TG5x$TotalTime = dataset$TG5x[ ,c(2:8,9)]
    new.data$TG5x$TotalTrials = dataset$TG5x[ ,c(2:8,10)]
    new.data$TG5x$Accuracy = dataset$TG5x[ ,c(2:8,11)]
    new.data$TG5x$Omission = dataset$TG5x[ ,c(2:8,12)]
    new.data$TG5x$Premature = dataset$TG5x[ ,c(2:8,13)]
    new.data$TG5x$Perseverative = dataset$TG5x[ ,c(2:8,14)]
    new.data$TG5x$RewardLat = dataset$TG5x[ ,c(2:8,15)]
    new.data$TG5x$CorrectLat = dataset$TG5x[ ,c(2:8,16)]
    
    new.data$TG3x$TotalTime = dataset$TG3x[ ,c(2:8,9)]
    new.data$TG3x$TotalTrials = dataset$TG3x[ ,c(2:8,10)]
    new.data$TG3x$Accuracy = dataset$TG3x[ ,c(2:8,11)]
    new.data$TG3x$Omission = dataset$TG3x[ ,c(2:8,12)]
    new.data$TG3x$Premature = dataset$TG3x[ ,c(2:8,13)]
    new.data$TG3x$Perseverative = dataset$TG3x[ ,c(2:8,14)]
    new.data$TG3x$RewardLat = dataset$TG3x[ ,c(2:8,15)]
    new.data$TG3x$CorrectLat = dataset$TG3x[ ,c(2:8,16)]
  }
  return(new.data)
}

# Data Format - Long to Wide ##
Data.Formatting.Function = function(dataset,m.value){
  for(a in 1:3){
    for(b in 1:m.value){
      temp.data = dataset[[a]][b]
      colnames(temp.data)[8] = 'Value'
      data.melt = melt(temp.data, id.vars = c("AnimalID","TestSite","Mouse.Strain","Genotype","Sex","Age.Months","Stimulus.Length"))
      data.melt$Age = as.character(data.melt$Age)
      data.melt$Value = as.numeric(data.melt$Value)
      data.cast = dcast(data.melt, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Stimulus.Length, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
      dataset[[a]][b] = data.cast
    }
  }
  return(dataset)
}

## Settings ##
options(scipen=50)
options(contrasts = c('contr.sum','contr.poly'))

## Read Data ##
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Pretrain QC Oct 12 2017 Updated.csv')
raw.data.acq = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Acquisition Aggregated QC Oct 12 2017 Updated.csv')
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

## Separate each Raw File by Strain ##
pretrain.separated.data = Strain.Separation.Function(raw.data.pretrain,0)
acq.separated.data = Strain.Separation.Function(raw.data.acq,0)
probe.separated.data = Strain.Separation.Function(raw.data.probe,1)

## Separate Probe Data by Measure / Get Specific Columns ##
pretrain.separated.measures = Measure.Separation.Function(pretrain.separated.data,0)
acq.separated.measures = Measure.Separation.Function(acq.separated.data,0)
probe.seperated.measures = Measure.Separation.Function(probe.separated.data,1)

## Format Data Long to Wide ##
pretrain.formatted.data = Data.Formatting.Function(pretrain.separated.measures,1)
acq.formatted.data = Data.Formatting.Function(acq.separated.measures,1)
probe.formatted.data = Data.Formatting.Function(probe.seperated.measures,8)

