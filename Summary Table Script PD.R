#################################################
#                                               #
#                                               #
#         Weston Anova Script                   #
#        Pairwise Discrimination                #
#                                               #
#                                               #
#################################################

## Library ##
library(plyr)
#library(tidyverse)
library(reshape2)
library(car)

## Functions List ##

# Separate Into Strains for Analysis #
Strain.Separation.Function = function(dataset,long.form=0){
  new.data = list()
  if(long.form == 0){
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), c(2:9)]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), c(2:9)]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="5xFAD"), c(2:9)] 
  }else if(long.form == 1){
    new.data$APP = dataset[which(dataset$Mouse.Strain=="APP-PS1"), c(2:11,15:18,49,81)]
    new.data$TG3x = dataset[which(dataset$Mouse.Strain=="3xTG-AD"), c(2:11,15:18,49,81)]
    new.data$TG5x = dataset[which(dataset$Mouse.Strain=="5xFAD"), c(2:11,15:18,49,81)]
  }
  return(new.data)
}

# Separate Into Separate Measures for Each File + Transform # #
Measure.Separation.Function = function(dataset, file.type=0){
  new.data = list()
  if(file.type == 1){
    new.data$APP$Sessions = dataset$APP[ ,c(1:8)]
    new.data$TG5x$Sessions = dataset$TG5x[ ,c(1:8)]
    new.data$TG3x$Sessions = dataset$TG3x[ ,c(1:8)]
  }else if(file.type == 2){
    new.data$APP$Sessions = dataset$APP[ ,c(1:8)]
    new.data$TG5x$Sessions = dataset$TG5x[ ,c(1:8)]
    new.data$TG3x$Sessions = dataset$TG3x[ ,c(1:8)]
  }else if(file.type == 3){
    new.data$APP$TotalTime = dataset$APP[ ,c(2:8,10,11)]
    new.data$APP$TotalTrials = dataset$APP[ ,c(2:8,10,12)]
    new.data$APP$Accuracy = dataset$APP[ ,c(2:8,10,14)]
    new.data$APP$Corrections = dataset$APP[ ,c(2:8,10,13)]
    new.data$APP$RewardLat = dataset$APP[ ,c(2:8,10,16)]
    new.data$APP$CorrectLat = dataset$APP[ ,c(2:8,10,15)]
    
    new.data$TG5x$TotalTime = dataset$TG5x[ ,c(2:8,10,11)]
    new.data$TG5x$TotalTrials = dataset$TG5x[ ,c(2:8,10,12)]
    new.data$TG5x$Accuracy = dataset$TG5x[ ,c(2:8,10,14)]
    new.data$TG5x$Corrections = dataset$TG5x[ ,c(2:8,10,13)]
    new.data$TG5x$RewardLat = dataset$TG5x[ ,c(2:8,10,16)]
    new.data$TG5x$CorrectLat = dataset$TG5x[ ,c(2:8,10,15)]
    
    new.data$TG3x$TotalTime = dataset$TG3x[ ,c(2:8,10,11)]
    new.data$TG3x$TotalTrials = dataset$TG3x[ ,c(2:8,10,12)]
    new.data$TG3x$Accuracy = dataset$TG3x[ ,c(2:8,10,14)]
    new.data$TG3x$Corrections = dataset$TG3x[ ,c(2:8,10,13)]
    new.data$TG3x$RewardLat = dataset$TG3x[ ,c(2:8,10,16)]
    new.data$TG3x$CorrectLat = dataset$TG3x[ ,c(2:8,10,15)]
  }
  return(new.data)
}

# Data Format - Long to Wide ##
Data.Formatting.Function = function(dataset,m.value,data.type){
  for(a in 1:3){
    for(b in 1:m.value){
      temp.data = as.data.frame(dataset[[a]][[b]])
      if(data.type == 3){
        colnames(temp.data)[9] = 'Value'
        temp.data = temp.data[which((temp.data$Task == 'PD Reversal 1') |(temp.data$Task == 'PD Reversal 2') | (temp.data$Task == 'PD Reversal 3')), ]
      }else{
        colnames(temp.data)[8] = 'Value' 
      }
      if(data.type == 1){
        data.cast = dcast(temp.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Task, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
      }else if(data.type == 2){
        data.cast = dcast(temp.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Task, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
      }else if(data.type == 3){
        data.cast = dcast(temp.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex ~ Age.Months + Task + Day, fun.aggregate = mean, na.rm=TRUE, value.var="Value")
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
  template.data = template.data[which((template.data$Task == 'PD Reversal 1') |(template.data$Task == 'PD Reversal 2') | (template.data$Task == 'PD Reversal 3')), ]
  idata = unique(template.data[c('Day')])
  idata = idata[order(idata), ]
  idata = as.factor(idata)
  return(idata)
}

# Calculate ANOVA Results #
Anova.Preparation.Main.Function = function(dataset,idata,age){
  for(a in 1:length(dataset)){
    for(b in 1:length(dataset[[a]])){
      data.file =dataset[[a]][[b]]
      data.file[ ,c(1,3)] = NULL
      if(age == 4){
        data.file = data.file[ ,c(1:3,4:13)]
        data.file = data.file[complete.cases(data.file), ]
        data.depend = data.file[ ,4:13]
        
      }else if(age == 7){
        data.file = data.file[ ,c(1:3,14:23)]
        data.file = data.file[complete.cases(data.file), ]
        data.depend = data.file[ ,4:13]
        
      }else if(age == 10){
        data.file = data.file[ ,c(1:3,24:33)]
        data.file = data.file[complete.cases(data.file), ]
        data.depend = data.file[ ,4:13]
      }
      data.lm = lm(as.matrix(data.depend) ~ 1 + Genotype * Sex, data=data.file)
      Day = as.data.frame(matrix(nrow=10,ncol=0))
      Day$Day = as.factor(idata)
      data.anova = Anova(data.lm, idata=Day,idesign=~Day, type="III")
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
      data.file = data.file[ ,c(4,5,9)]
      data.file = data.file[complete.cases(data.file), ]
      data.depend = data.file[ ,3]
      data.lm = lm(as.matrix(data.depend) ~ 1+ Genotype * Sex, data=data.file)
      data.anova = Anova(data.lm, type="III")
      dataset[[a]][[b]] = data.anova
    }
  }
  return(dataset)
}

Anova.Preparation.Acquisition.Function = function(dataset,age){
  final.dataset = list()
  for(a in 1:length(dataset)){
    for(b in 1:length(dataset[[a]])){
      data.file =dataset[[a]][[b]]
      data.file = data.file[complete.cases(data.file), ]
      data.file[ ,c(1:3)] = NULL
      if(age == 4){
        data.depend = data.file[ ,3]
      }else if(age == 7){
        data.depend = data.file[ ,4]
      }else if(age == 10){
        data.depend = data.file[ ,5]
      }
      data.lm = lm(as.matrix(data.depend) ~ 1+ Genotype * Sex, data=data.file)
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
raw.data.pretrain = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PD\\Weston PD Pretrain QC Mar 26 2018.csv')
raw.data.acq = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PD\\Weston PD Acquisition Aggregated Mar 26 2018.csv')
raw.data.main = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PD\\Weston PD Basline Reversal QC Mar 26 2018.csv')
raw.data.acq[ ,2] = NULL
raw.data.main$No.Correction.Trials[is.na(raw.data.main$No.Correction.Trials)] = 0
## Separate each Raw File by Strain ##
pretrain.separated.data = Strain.Separation.Function(raw.data.pretrain,0)
acq.separated.data = Strain.Separation.Function(raw.data.acq,0)
main.separated.data = Strain.Separation.Function(raw.data.main,1)


## Separate Probe Data by Measure / Get Specific Columns ##
pretrain.separated.measures = Measure.Separation.Function(pretrain.separated.data,1)
acq.separated.measures = Measure.Separation.Function(acq.separated.data,2)
main.separated.measures = Measure.Separation.Function(main.separated.data,3)

## Remove Baseline Protocol From Main Task ##

## Format Data Long to Wide ##
pretrain.formatted.data = Data.Formatting.Function(pretrain.separated.measures,1,1)
acq.formatted.data = Data.Formatting.Function(acq.separated.measures,1,1)
main.formatted.data = Data.Formatting.Function(main.separated.measures,6,3)

## Gather iData for Repeated Measures ##
main.idata = iData.Generate.Function(main.separated.measures)

## Conduct ANOVA ##
pretrain.anova = Anova.Preparation.Pretrain.Function(pretrain.formatted.data)
acq.4month.anova = Anova.Preparation.Acquisition.Function(acq.formatted.data,4)
acq.7month.anova = Anova.Preparation.Acquisition.Function(acq.formatted.data,7)
acq.10month.anova = Anova.Preparation.Acquisition.Function(acq.formatted.data,10)
main.4month.anova = Anova.Preparation.Main.Function(main.formatted.data,main.idata,4)
main.7month.anova = Anova.Preparation.Main.Function(main.formatted.data,main.idata,7)
main.10month.anova = Anova.Preparation.Main.Function(main.formatted.data, main.idata, 10)

## Probe Anova Prep ##
probe.anova = list()
probe.anova$APP$Month4 = main.4month.anova$APP
probe.anova$APP$Month7 = main.7month.anova$APP
probe.anova$APP$Month10 = main.10month.anova$APP

probe.anova$TG5x$Month4 = main.4month.anova$TG5x
probe.anova$TG5x$Month7 = main.7month.anova$TG5x
probe.anova$TG5x$Month10 = main.10month.anova$TG5x

probe.anova$TG3x$Month4 = main.4month.anova$TG3x
probe.anova$TG3x$Month7 = main.7month.anova$TG3x
probe.anova$TG3x$Month10 = main.10month.anova$TG3x

## Prepare ANOVA Table ##
strain.list = as.vector(names(probe.anova))
measure.list = as.vector(names(probe.anova$APP$Month4))
measure.list = measure.list[c(3:6)]
age.list = as.vector(names(probe.anova$APP))
template.file = summary(probe.anova[[1]][[1]][[1]], multivariate=FALSE)
template.rownames = rownames(template.file[[4]])
template.rownames = template.rownames[2:length(template.rownames)]
template.pvalads = rownames(template.file[[5]])

hm.rownames = template.rownames
hm.rownames = gsub('Genotype','G',hm.rownames)
hm.rownames = gsub('Sex','Sx',hm.rownames)
hm.rownames = gsub('Day','D',hm.rownames)
hm.rownames = gsub(':','*',hm.rownames)

strain.count = length(strain.list)
age.count = length(age.list)
measure.count = length(measure.list)
analysis.count = length(hm.rownames)

map.list = list()
summary.table = as.data.frame(matrix(nrow=(strain.count*age.count*measure.count*analysis.count),ncol=9))
colnames(summary.table) = c('Strain','Age','Measure','Analysis','df1','df2','F','p','partial eta^2')

for(a in 1:length(strain.list)){
  for(b in 1:length(age.list)){
    for(d in 1:length(measure.list)){
      temp.summary = summary(probe.anova[[strain.list[a]]][[age.list[b]]][[measure.list[d]]], multivariate=FALSE)
      temp.main = temp.summary[[4]]
      temp.pvalad = temp.summary[[5]]
      for(c in 1:length(template.pvalads)){
        temp.main[which(rownames(temp.main) == template.pvalads[c]),6] = temp.pvalad[which(rownames(temp.pvalad) == template.pvalads[c]),2]
      }
      for(c in 1:length(hm.rownames)){
        row.modifier = (((a-1)*age.count*measure.count*analysis.count) + ((b-1)*measure.count*analysis.count) + ((d-1)*analysis.count) + c)
        partial.eta = temp.main[(c+1),1] / (temp.main[(c+1),1] + temp.main[(c+1),3])
        summary.table[row.modifier,1] = strain.list[a]
        summary.table[row.modifier,2] = age.list[b]
        summary.table[row.modifier,3] = measure.list[d]
        summary.table[row.modifier,4] = hm.rownames[c]
        summary.table[row.modifier,5] = temp.main[(c+1),2]
        summary.table[row.modifier,6] = temp.main[(c+1),4]
        summary.table[row.modifier,7] = round(temp.main[(c+1),5],digits=2)
        summary.table[row.modifier,8] = round(temp.main[(c+1),6],digits=3)
        summary.table[row.modifier,9] = round(partial.eta, digits = 2)
      }
    }
  }
}
map.list[['Summary']] = summary.table