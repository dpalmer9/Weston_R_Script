#################################################
#                                               #
#                                               #
#         Weston Anova Script                   #
#             PAL                               #
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
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), c(2:10)]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="5xFAD"), c(2:10)]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), c(2:10)] 
  }else if(long.form == 1){
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), c(2:10,13:16,239,315)]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="5xFAD"), c(2:10,13:16,239,315)]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), c(2:10,13:16,239,315)]
  }
  return(new.data)
}

# Separate Into Separate Measures for Each File + Transform # #
Measure.Separation.Function = function(dataset, long.form=0){
  new.data = list()
  if(long.form == 0){
    new.data$APP$Sessions = dataset$APP
    new.data$TG5x$Sessions = dataset$TG5x
    new.data$TG3x$Sessions = dataset$TG3x
  }else if(long.form == 1){
    new.data$APP$TotalTime = dataset$APP[ ,c(1:8,10)]
    new.data$APP$TotalTrials = dataset$APP[ ,c(1:8,11)]
    new.data$APP$Accuracy = dataset$APP[ ,c(1:8,13)]
    new.data$APP$Corrections = dataset$APP[ ,c(1:8,12)]
    new.data$APP$RewardLat = dataset$APP[ ,c(1:8,15)]
    new.data$APP$CorrectLat = dataset$APP[ ,c(1:8,14)]
    
    new.data$TG5x$TotalTime = dataset$TG5x[ ,c(1:8,10)]
    new.data$TG5x$TotalTrials = dataset$TG5x[ ,c(1:8,11)]
    new.data$TG5x$Accuracy = dataset$TG5x[ ,c(1:8,13)]
    new.data$TG5x$Corrections = dataset$TG5x[ ,c(1:8,12)]
    new.data$TG5x$RewardLat = dataset$TG5x[ ,c(1:8,15)]
    new.data$TG5x$CorrectLat = dataset$TG5x[ ,c(1:8,14)]
    
    new.data$TG3x$TotalTime = dataset$TG3x[ ,c(1:8,10)]
    new.data$TG3x$TotalTrials = dataset$TG3x[ ,c(1:8,11)]
    new.data$TG3x$Accuracy = dataset$TG3x[ ,c(1:8,13)]
    new.data$TG3x$Corrections = dataset$TG3x[ ,c(1:8,12)]
    new.data$TG3x$RewardLat = dataset$TG3x[ ,c(1:8,15)]
    new.data$TG3x$CorrectLat = dataset$TG3x[ ,c(1:8,14)]
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

raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PAL\\Weston PAL Pretraining Oct 12 2017 Updated.csv')
raw.data.acq = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PAL\\Weston PAL Acquisition Aggregated Oct 12 2017 Updated.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PAL\\Weston PAL Main Task Aggregated Oct 12 2017 Updated.csv')

colnames(raw.data.pretrain)[6] = "Sex"
## Separate each Raw File by Strain ##
pretrain.separated.data = Strain.Separation.Function(raw.data.pretrain,0)
acq.separated.data = Strain.Separation.Function(raw.data.acq,0)
main.separated.data = Strain.Separation.Function(raw.data.main,1)


## Separate Probe Data by Measure / Get Specific Columns ##
pretrain.separated.measures = Measure.Separation.Function(pretrain.separated.data,0)
acq.separated.measures = Measure.Separation.Function(acq.separated.data,0)
main.separated.measures = Measure.Separation.Function(main.separated.data,1)

## Format Data Long to Wide ##
pretrain.formatted.data = Data.Formatting.Function(pretrain.separated.measures,1,0)
acq.formatted.data = Data.Formatting.Function(acq.separated.measures,1,1)
main.formatted.data = Data.Formatting.Function(main.separated.measures,8,1)
