## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(mice)
#library(lme4)
#library(nlm3)


## Parameters ##


## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceProbe.csv')
#raw.data = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\FiveChoiceProbe.csv')

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
data.3x.list.probe = mget(ls(pattern="data.3x.probe."))

data.5x.probe.totaltime = data.5x.probe[c(1:7,8)]
data.5x.probe.totaltrial = data.5x.probe[c(1:7,9)]
data.5x.probe.acc = data.5x.probe[c(1:7,10)]
data.5x.probe.omission = data.5x.probe[c(1:7,11)]
data.5x.probe.premature = data.5x.probe[c(1:7,12)]
data.5x.probe.persev = data.5x.probe[c(1:7,13)]
data.5x.probe.rewardlat = data.5x.probe[c(1:7,14)]
data.5x.probe.corrlat = data.5x.probe[c(1:7,15)]
data.5x.list.probe = mget(ls(pattern="data.5x.probe."))

data.app.probe.totaltime = data.app.probe[c(1:7,8)]
data.app.probe.totaltrial = data.app.probe[c(1:7,9)]
data.app.probe.acc = data.app.probe[c(1:7,10)]
data.app.probe.omission = data.app.probe[c(1:7,11)]
data.app.probe.premature = data.app.probe[c(1:7,12)]
data.app.probe.persev = data.app.probe[c(1:7,13)]
data.app.probe.rewardlat = data.app.probe[c(1:7,14)]
data.app.probe.corrlat = data.app.probe[c(1:7,15)]
data.app.list.probe = mget(ls(pattern="data.app.probe."))

## Format Data Long to Wide ##
Data.Spread.Function <- function(dataset){
  data.melt = melt(dataset, id.vars = c("AnimalID","Site","Strain","Genotype","Gender","Age","ProbeDuration"))
  data.melt$Age = as.character(data.melt$Age)
  data.melt$Age[data.melt$Age == "4"] = "04"
  data.melt$Age[data.melt$Age == "7"] = "07"
  data.melt$Age[data.melt$Age == "13_15"] = "13"
  data.melt$value = as.numeric(data.melt$value)
  data.cast = dcast(data.melt, AnimalID + Site + Strain + Genotype + Gender ~ variable + Age + ProbeDuration, fun.aggregate = mean, na.rm=TRUE, value.var="value")
  return(data.cast)
}

data.3x.list.melt = lapply(data.3x.list.probe, Data.Spread.Function)
for(a in 1:length(data.3x.list.probe)){
  assign(names(data.3x.list.probe)[a], as.data.frame(data.3x.list.probe[a]))
}
data.5x.list.melt = lapply(data.5x.list.probe, Data.Spread.Function)
for(a in 1:length(data.5x.list.probe)){
  assign(names(data.5x.list.probe)[a], as.data.frame(data.5x.list.probe[a]))
}
data.app.list.melt = lapply(data.app.list.probe, Data.Spread.Function)
for(a in 1:length(data.app.list.probe)){
  assign(names(data.app.list.probe)[a], as.data.frame(data.app.list.probe[a]))
}
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
  data.lm = lm(as.matrix(data.depend) ~ 1+ Site * Genotype * Gender, data=dataset)
  data.anova = Anova(data.lm, idata=idata,idesign=~Age*ProbeDuration, type="III")
  return(data.anova)
}

## Generate lm for each file - REML Function ##
Data.GenerateMI.Function <- function(dataset){
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Gender')
  dataset$Strain = NULL
  dataset$AnimalID = NULL
  dataset.mi = mice(dataset, m=5, maxit=50, method='pmm',seed=500)
  dataset.1 = complete(dataset.mi, 1)
  data.lm.1 = lm(as.matrix(dataset[4:length(colnames(dataset.1))]) ~ 1 + Site * Genotype * Gender, data=dataset.1)
  dataset.2 = complete(dataset.mi, 2)
  data.lm.2 = lm(as.matrix(dataset[4:length(colnames(dataset.2))]) ~ 1 + Site * Genotype * Gender, data=dataset.2)
  dataset.3 = complete(dataset.mi, 3)
  data.lm.3 = lm(as.matrix(dataset[4:length(colnames(dataset.3))]) ~ 1 + Site * Genotype * Gender, data=dataset.3)
  dataset.4 = complete(dataset.mi, 4)
  data.lm.4 = lm(as.matrix(dataset[4:length(colnames(dataset.4))]) ~ 1 + Site * Genotype * Gender, data=dataset.4)
  dataset.5 = complete(dataset.mi, 5)
  data.lm.5 = lm(as.matrix(dataset[4:length(colnames(dataset.5))]) ~ 1 + Site * Genotype * Gender, data=dataset.5)
  data.mira = as.mira(list(data.lm.1,data.lm.2,data.lm.3,data.lm.4,data.lm.5))
  dataset.pool = pool(data.mira)
  dataset.fit = with(data= dataset.mi, exp = lm(as.matrix(dataset[4:length(colnames(dataset.mi))]) ~ 1 + Site * Genotype * Gender))
  dataset.pool2 = pool(dataset.fit)
  testlist = list(data.lm.1,data.lm.2,data.lm.3,data.lm.4,data.lm.5)
  dataset.anova = miceadds::mi.anova(mi.res=dataset.mi, formula = "as.matrix(dataset[4:length(colnames(dataset.mi))]) ~ 1 + Site * Genotype * Gender", type=3)
  return(data.lm.5)
}

Data.GenerateMI.Function2 <- function(dataset){
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Gender')
  dataset.melt = melt(dataset, id.vars = c('AnimalID','Site','Strain','Genotype','Gender'))
  measure.list = as.data.frame(stringr::str_split_fixed(dataset.melt$variable, "_", 3))
  dataset.melt$Age = measure.list$V2
  dataset.melt$Age = as.factor(dataset.melt$Age)
  dataset.melt$ProbeDur = measure.list$V3
  dataset.melt$ProbeDur = as.factor(dataset.melt$ProbeDur)
  dataset.fixed = dataset.melt[ ,c(2,4,5,8,9,7)]
  colnames(dataset.fixed) = c('Site','Genotype','Sex','Age','ProbeDur','Measure')
  dataset.mi = mice(dataset.fixed, m=5, maxit=50, method='pmm',seed=500)
  #dataset.fit = with(data= dataset.mi, exp = lmer(Measure ~ Genotype + Sex + Site + (1|Age) + (1|ProbeDur)))
  #dataset.pool = pool(dataset.fit)
  return(dataset.mi)
}

## Pass ANOVA function through lists ##
data.3x.list.lmanova = lapply(data.3x.list.melt, Data.GenerateLM.Function, idata=data.3x.idata)

data.5x.list.lmanova = lapply(data.5x.list.melt, Data.GenerateLM.Function, idata=data.5x.idata)

data.app.list.lmanova = lapply(data.app.list.melt, Data.GenerateLM.Function, idata=data.app.idata)

## Vigilance Configuration ##
Data.SeparateVigilance.Function <- function(dataset, genotype, age="NA", sex="NA"){
  dataset$Age..months. = as.character(dataset$Age..months.)
  if(genotype == "w"){
    new.dataset = dataset[which(dataset$Genotype == "w"), ]
  }else if(genotype == "t"){
    new.dataset = dataset[which(dataset$Genotype == "t"), ]
  }
  if(age == "4"){
    new.dataset = dataset[which(new.dataset$Age..months. == "4"), ]
  }else if(age == "7"){
    new.dataset = dataset[which(new.dataset$Age..months. == "7"), ]
  }else if(age == "10"){
    new.dataset = dataset[which(new.dataset$Age..months. == "10"), ]
  }
  if(sex == "F"){
    new.dataset = dataset[which(new.dataset$Gender == "F"), ]
  }else if(sex == "M"){
    new.dataset = dataset[which(new.dataset$Gender == "M"), ]
  }
  datalist = list()
  datalist$acc = new.dataset[ ,c(3,9,29,39,49,59,69)]
  colnames(datalist$acc) = c('AnimalID','ProbeDur','Bin10','Bin20','Bin30','Bin40','Bin50')
  datalist$acc = melt(datalist$acc, id.vars = c('AnimalID','ProbeDur'))
  datalist$acc = dcast(datalist$acc, AnimalID ~ ProbeDur + variable, fun.aggregate = mean, na.rm=TRUE)
  datalist$omit = new.dataset[ ,c(3,9,81,91,101,111,121)]
  colnames(datalist$omit) = c('AnimalID','ProbeDur','Bin10','Bin20','Bin30','Bin40','Bin50')
  datalist$omit = melt(datalist$omit, id.vars = c('AnimalID','ProbeDur'))
  datalist$omit = dcast(datalist$omit, AnimalID ~ ProbeDur + variable, fun.aggregate = mean, na.rm=TRUE)
  return(datalist)
  
  #datalist = 
} 

vig.3x.wt.m.4 = Data.SeparateVigilance.Function(data.3x, "w", "4","M")
vig.3x.wt.m.7 = Data.SeparateVigilance.Function(data.3x, "w", "7","M")
vig.3x.wt.m.10 = Data.SeparateVigilance.Function(data.3x, "w", "10","M")

vig.3x.tg.m.4 = Data.SeparateVigilance.Function(data.3x, "t", "4","M")
vig.3x.tg.m.7 = Data.SeparateVigilance.Function(data.3x, "t", "7","M")
vig.3x.tg.m.10 = Data.SeparateVigilance.Function(data.3x, "t", "10","M")

vig.5x.wt.m.4 = Data.SeparateVigilance.Function(data.5x, "w", "4","M")
vig.5x.wt.m.7 = Data.SeparateVigilance.Function(data.5x, "w", "7","M")
vig.5x.wt.m.10 = Data.SeparateVigilance.Function(data.5x, "w", "10","M")

vig.5x.tg.m.4 = Data.SeparateVigilance.Function(data.5x, "t", "4","M")
vig.5x.tg.m.7 = Data.SeparateVigilance.Function(data.5x, "t", "7","M")
vig.5x.tg.m.10 = Data.SeparateVigilance.Function(data.5x, "t", "10","M")

vig.app.wt.m.4 = Data.SeparateVigilance.Function(data.app, "w", "4","M")
vig.app.wt.m.7 = Data.SeparateVigilance.Function(data.app, "w", "7","M")
vig.app.wt.m.10 = Data.SeparateVigilance.Function(data.app, "w", "10","M")

vig.app.tg.m.4 = Data.SeparateVigilance.Function(data.app, "t", "4","M")
vig.app.tg.m.7 = Data.SeparateVigilance.Function(data.app, "t", "7","M")
vig.app.tg.m.10 = Data.SeparateVigilance.Function(data.app, "t", "10","M")

vig.3x.wt.f.4 = Data.SeparateVigilance.Function(data.3x, "w", "4","F")
vig.3x.wt.f.7 = Data.SeparateVigilance.Function(data.3x, "w", "7","F")
vig.3x.wt.f.10 = Data.SeparateVigilance.Function(data.3x, "w", "10","F")

vig.3x.tg.f.4 = Data.SeparateVigilance.Function(data.3x, "t", "4","F")
vig.3x.tg.f.7 = Data.SeparateVigilance.Function(data.3x, "t", "7","F")
vig.3x.tg.f.10 = Data.SeparateVigilance.Function(data.3x, "t", "10","F")

vig.5x.wt.f.4 = Data.SeparateVigilance.Function(data.5x, "w", "4","F")
vig.5x.wt.f.7 = Data.SeparateVigilance.Function(data.5x, "w", "7","F")
vig.5x.wt.f.10 = Data.SeparateVigilance.Function(data.5x, "w", "10","F")

vig.5x.tg.f.4 = Data.SeparateVigilance.Function(data.5x, "t", "4","F")
vig.5x.tg.f.7 = Data.SeparateVigilance.Function(data.5x, "t", "7","F")
vig.5x.tg.f.10 = Data.SeparateVigilance.Function(data.5x, "t", "10","F")

vig.app.wt.f.4 = Data.SeparateVigilance.Function(data.app, "w", "4","F")
vig.app.wt.f.7 = Data.SeparateVigilance.Function(data.app, "w", "7","F")
vig.app.wt.f.10 = Data.SeparateVigilance.Function(data.app, "w", "10","F")

vig.app.tg.f.4 = Data.SeparateVigilance.Function(data.app, "t", "4","F")
vig.app.tg.f.7 = Data.SeparateVigilance.Function(data.app, "t", "7","F")
vig.app.tg.f.10 = Data.SeparateVigilance.Function(data.app, "t", "10","F")


