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
Data.Formatting.Function = function(dataset,m.value,datatype){
  for(a in 1:3){
    for(b in 1:m.value){
      temp.data = as.data.frame(dataset[[a]][[b]])
      colnames(temp.data)[8] = 'Value'
      if(isTRUE(datatype == 0)){
        data.cast = dcast(temp.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Task, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
      }else if(isTRUE(datatype == 1)){
        data.cast = dcast(temp.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Stimulus.Length, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
      }
      for(c in 6:ncol(data.cast)){
        colnames(data.cast)[c] = paste('Data',colnames(data.cast)[c],sep=".")
      }
      dataset[[a]][[b]] = as.data.frame(data.cast)
    }
  }
  return(dataset)
}

# Generate iData for Repeated Measure Design (Probe Only) #
iData.Generate.Function = function(dataset){
  template.data = dataset[[1]][[1]]
  idata = unique(template.data[c('Age.Months','Stimulus.Length')])
  idata = idata[order(idata$Age.Months,idata$Stimulus.Length), ]
  idata$Age.Months = as.factor(idata$Age.Months)
  idata$Stimulus.Length = as.factor(idata$Stimulus.Length)
  return(idata)
}

# Calculate ANOVA Results #
Anova.Preparation.Probe.Function = function(dataset,idata){
  for(a in 1:length(dataset)){
    for(b in 1:length(dataset[[a]])){
      data.file =dataset[[a]][[b]]
      data.file = data.file[complete.cases(data.file), ]
      data.file[ ,c(1,3)] = NULL
      data.depend = data.file[ ,4:ncol(data.file)]
      data.lm = lm(as.matrix(data.depend) ~ 1+ TestSite * Genotype * Sex, data=data.file)
      data.anova = Anova(data.lm, idata=idata,idesign=~Age.Months*Stimulus.Length, type="III")
      dataset[[a]][[b]] = data.anova
    }
  }
  return(dataset)
}

Anova.Preparation.Pretrain.Function = function(dataset){
  final.dataset = list()
  for(a in 1:length(dataset)){
    for(b in 1:length(dataset[[a]])){
      data.file =dataset[[a]][[b]]
      data.file = data.file[complete.cases(data.file), ]
      data.file[ ,c(1,3)] = NULL
      data.depend = data.file[ ,7]
      data.lm = lm(as.matrix(data.depend) ~ 1+ TestSite * Genotype * Sex, data=data.file)
      data.anova = Anova(data.lm, type="III")
      dataset[[a]][[b]] = data.anova
    }
  }
  return(dataset)
}

Anova.Preparation.Acquisition.Function = function(dataset,trainingstim){
  final.dataset = list()
  for(a in 1:length(dataset)){
    for(b in 1:length(dataset[[a]])){
      data.file =dataset[[a]][[b]]
      data.file = data.file[complete.cases(data.file), ]
      data.file[ ,c(1,3)] = NULL
      if(trainingstim == 4){
        data.depend = data.file[ ,5]
      }else if(trainingstim == 2){
        data.depend = data.file[ ,4]
      }
      data.lm = lm(as.matrix(data.depend) ~ 1+ TestSite * Genotype * Sex, data=data.file)
      data.anova = Anova(data.lm, type="III")
      dataset[[a]][[b]] = data.anova
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

raw.data.pretrain = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Pretrain QC Oct 12 2017 Updated.csv')
raw.data.acq = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Acquisition Aggregated QC Oct 12 2017 Updated.csv')
raw.data.probe = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Pretrain QC Oct 12 2017 Updated.csv')
raw.data.acq = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Acquisition Aggregated QC Oct 12 2017 Updated.csv')
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

## Separate each Raw File by Strain ##
pretrain.separated.data = Strain.Separation.Function(raw.data.pretrain,0)
acq.separated.data = Strain.Separation.Function(raw.data.acq,0)
probe.separated.data = Strain.Separation.Function(raw.data.probe,1)

## Separate Probe Data by Measure / Get Specific Columns ##
pretrain.separated.measures = Measure.Separation.Function(pretrain.separated.data,0)
acq.separated.measures = Measure.Separation.Function(acq.separated.data,0)
probe.separated.measures = Measure.Separation.Function(probe.separated.data,1)

## Format Data Long to Wide ##
pretrain.formatted.data = Data.Formatting.Function(pretrain.separated.measures,1,0)
acq.formatted.data = Data.Formatting.Function(acq.separated.measures,1,1)
probe.formatted.data = Data.Formatting.Function(probe.separated.measures,8,1)

## Gather iData for Repeated Measures ##
probe.idata = iData.Generate.Function(probe.separated.measures)

## Conduct ANOVA ##
pretrain.anova = Anova.Preparation.Pretrain.Function(pretrain.formatted.data)
acq.4.anova = Anova.Preparation.Acquisition.Function(acq.formatted.data,4)
acq.2.anova = Anova.Preparation.Acquisition.Function(acq.formatted.data,2)
probe.anova = Anova.Preparation.Probe.Function(probe.formatted.data,probe.idata)
