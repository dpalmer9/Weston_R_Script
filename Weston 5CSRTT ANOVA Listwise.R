## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
library(nlm3)


## Parameters ##


## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceProbe.csv')


## Separate By Strain ##
data.3x = raw.data[which(raw.data$Mouse.strain=="3XTG"), ]
data.5x = raw.data[which(raw.data$Mouse.strain=="5FAD"), ]
data.app = raw.data[which(raw.data$Mouse.strain=="APPPS1"), ]

## Partion Data to Probe ##

data.3x.probe = data.3x[c(3:9,14:19,176,228)]
data.5x.probe = data.5x[c(3:9,14:19,176,228)]
data.app.probe = data.app[c(3:9,14:19,176,228)]

## Fix Column Names ##
Data.ColumnName.Fix <- function(dataset){
  colnames(dataset) <- c('AnimalID','Site','Strain','Genotype','Gender','Age','ProbeDuration','TotalTime','TotalTrial','Accuracy','Omissions','PrematureResponse','PerseverativeResponse','RewardLatency','CorrectLatency')
  return(dataset)
}

data.3x.probe = Data.ColumnName.Fix(data.3x.probe)
data.5x.probe = Data.ColumnName.Fix(data.5x.probe)
data.app.probe = Data.ColumnName.Fix(data.app.probe)

## Separate Data Into Var ##
data.3x.probe.totaltime = data.3x.probe[c(1:7,8)]
data.3x.probe.totaltrial = data.3x.probe[c(1:7,9)]
data.3x.probe.acc = data.3x.probe[c(1:7,10)]
data.3x.probe.omission = data.3x.probe[c(1:7,11)]
data.3x.probe.premature = data.3x.probe[c(1:7,12)]
data.3x.probe.persev = data.3x.probe[c(1:7,13)]
data.3x.probe.rewardlat = data.3x.probe[c(1:7,14)]
data.3x.probe.corrlat = data.3x.probe[c(1:7,15)]
data.3x.list.probe = ls(pattern="data.3x.probe.")

data.5x.probe.totaltime = data.5x.probe[c(1:7,8)]
data.5x.probe.totaltrial = data.5x.probe[c(1:7,9)]
data.5x.probe.acc = data.5x.probe[c(1:7,10)]
data.5x.probe.omission = data.5x.probe[c(1:7,11)]
data.5x.probe.premature = data.5x.probe[c(1:7,12)]
data.5x.probe.persev = data.5x.probe[c(1:7,13)]
data.5x.probe.rewardlat = data.5x.probe[c(1:7,14)]
data.5x.probe.corrlat = data.5x.probe[c(1:7,15)]
data.5x.list.probe = ls(pattern="data.5x.probe.")

data.app.probe.totaltime = data.app.probe[c(1:7,8)]
data.app.probe.totaltrial = data.app.probe[c(1:7,9)]
data.app.probe.acc = data.app.probe[c(1:7,10)]
data.app.probe.omission = data.app.probe[c(1:7,11)]
data.app.probe.premature = data.app.probe[c(1:7,12)]
data.app.probe.persev = data.app.probe[c(1:7,13)]
data.app.probe.rewardlat = data.app.probe[c(1:7,14)]
data.app.probe.corrlat = data.app.probe[c(1:7,15)]
data.app.list.probe = ls(pattern="data.app.probe.")

## Format Data Long to Wide ##
Data.SpreadLM.Function <- function(dataset){
  data.melt = melt(dataset, id.vars = c("AnimalID","Site","Strain","Genotype","Gender","Age","ProbeDuration"))
  data.melt$Age = as.character(data.melt$Age)
  data.melt$Age[data.melt$Age == "4"] = "04"
  data.melt$Age[data.melt$Age == "7"] = "07"
  data.melt$Age[data.melt$Age == "13_15"] = "13"
  data.melt$value = as.numeric(data.melt$value)
  data.cast = dcast(data.melt, AnimalID + Site + Strain + Genotype + Gender ~ variable + Age + ProbeDuration, fun.aggregate = mean, na.rm=TRUE, value.var="value")
  return(data.cast)
}

Data.SpreadREML.Function <- function(dataset){
  data.melt = melt(dataset, id.vars = c("AnimalID","Site","Strain","Genotype","Gender","Age","ProbeDuration"), na.rm=FALSE)
  data.melt$Age = as.character(data.melt$Age)
  data.melt$Age[data.melt$Age == "4"] = "04"
  data.melt$Age[data.melt$Age == "7"] = "07"
  data.melt$Age[data.melt$Age == "13_15"] = "13"
  data.melt$value = as.numeric(data.melt$value)
  return(data.melt)
}

data.3x.list.probe.lm = lapply(mget(data.3x.list.probe), Data.SpreadLM.Function)
#for(a in 1:length(data.3x.list.probe)){
  #assign(names(data.3x.list.probe)[a], as.data.frame(data.3x.list.probe[a]))
#}
data.5x.list.probe.lm = lapply(mget(data.5x.list.probe), Data.SpreadLM.Function)
#for(a in 1:length(data.5x.list.probe)){
  #assign(names(data.5x.list.probe)[a], as.data.frame(data.5x.list.probe[a]))
#}
data.app.list.probe.lm = lapply(mget(data.app.list.probe), Data.SpreadLM.Function)
#for(a in 1:length(data.app.list.probe)){
  #assign(names(data.app.list.probe)[a], as.data.frame(data.app.list.probe[a]))
#}

data.3x.list.probe.reml = lapply(mget(data.3x.list.probe), Data.SpreadREML.Function)
data.5x.list.probe.reml = lapply(mget(data.5x.list.probe), Data.SpreadREML.Function)
data.app.list.probe.reml = lapply(mget(data.app.list.probe), Data.SpreadREML.Function)
## Generate iData ##
Data.CreateiData.Function <- function(dataset){
  idata = unique(dataset[c('Age','ProbeDuration')])
  idata$Age = as.character(idata$Age)
  idata$Age[idata$Age == "4"] = "04"
  idata$Age[idata$Age == "7"] = "07"
  idata$Age[idata$Age == "13_15"] = "13"
  idata = idata[order(idata$Age,idata$ProbeDuration), ]
  return(idata)
}

data.3x.idata = Data.CreateiData.Function(data.3x.probe)
data.3x.idata$Age = as.factor(data.3x.idata$Age)
data.3x.idata$ProbeDuration = as.factor(data.3x.idata$ProbeDuration)

data.5x.idata = Data.CreateiData.Function(data.5x.probe)
data.5x.idata$Age = as.factor(data.5x.idata$Age)
data.5x.idata$ProbeDuration = as.factor(data.5x.idata$ProbeDuration)

data.app.idata = Data.CreateiData.Function(data.app.probe)
data.app.idata$Age = as.factor(data.app.idata$Age)
data.app.idata$ProbeDuration = as.factor(data.app.idata$ProbeDuration)


## Generate lm for each file - Listwise Function ##
Data.GenerateLM.Function <- function(dataset,idata){
  dataset = dataset[complete.cases(dataset), ]
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Gender')
  dataset$Strain = NULL
  dataset$AnimalID = NULL
  data.depend = dataset[4:length(colnames(dataset))]
  data.lm = lm(as.matrix(data.depend) ~ 1 * Site * Genotype * Gender, data=dataset)
  data.anova = Anova(data.lm, idata=idata,idesign=~Age*ProbeDuration, type="III")
  return(data.anova)
}
 
## Generate lm for each file - REML Function ##
Data.GenerateREML.Function <- function(dataset){
  dataset$variable = NULL
  dataset$Strain = NULL
  data.reml = lmer(value~Site + Genotype + Gender + (1|Age) + (1|ProbeDuration), data=dataset, REML = FALSE)
  data.anova = Anova(data.reml)
  return(data.anova)
}

## Pass ANOVA functions through lists ##
data.3x.list.probe.lm = lapply(data.3x.list.probe.lm, Data.GenerateLM.Function, idata=data.3x.idata)

data.5x.list.probe.lm = lapply(data.5x.list.probe.lm, Data.GenerateLM.Function, idata=data.5x.idata)

data.app.list.probe.lm = lapply(data.app.list.probe.lm, Data.GenerateLM.Function, idata=data.app.idata)
