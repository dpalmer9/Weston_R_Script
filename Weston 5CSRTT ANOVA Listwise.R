## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(mice)
library(lme4)
library(nlme)


## Parameters ##


## Acquire data ##
options(scipen=50)
options(contrasts = c('contr.sum','contr.poly'))
#raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceProbe.csv')
raw.data = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\FiveChoiceProbe.csv')

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
data.app.probe = data.app.probe[which(data.app.probe$Age != "13_15"), ]

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
  data.cast = dcast(data.melt, AnimalID + Site + Strain + Genotype + Gender ~ Age + ProbeDuration, fun.aggregate = mean, na.rm=TRUE, value.var="value")
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
Data.GenerateMI.Function <- function(dataset,idata){
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Gender')
  dataset$Strain = NULL
  dataset$AnimalID = NULL
  dataset.mi = mice(dataset, m=5, maxit=50, method='pmm',seed=500)
  dataset.1 = complete(dataset.mi, 1)
  data.lm.1 = lm(as.matrix(dataset.1[4:length(colnames(dataset.1))]) ~ 1 + Site * Genotype * Gender, data=dataset.1)
  dataset.anova.1 = Anova(data.lm.1, idata=idata,idesign=~Age*ProbeDuration, type="III")
  dataset.summary.1 = summary(dataset.anova.1, multivariate=FALSE)
  dataset.anovasheet.1 = as.matrix(dataset.summary.1$univariate.tests)
  
  dataset.2 = complete(dataset.mi, 2)
  data.lm.2 = lm(as.matrix(dataset.2[4:length(colnames(dataset.2))]) ~ 1 + Site * Genotype * Gender, data=dataset.2)
  dataset.anova.2 = Anova(data.lm.2, idata=idata,idesign=~Age*ProbeDuration, type="III")
  dataset.summary.2 = summary(dataset.anova.2, multivariate=FALSE)
  dataset.anovasheet.2 = as.matrix(dataset.summary.2$univariate.tests)
  
  dataset.3 = complete(dataset.mi, 3)
  data.lm.3 = lm(as.matrix(dataset.3[4:length(colnames(dataset.3))]) ~ 1 + Site * Genotype * Gender, data=dataset.3)
  dataset.anova.3 = Anova(data.lm.3, idata=idata,idesign=~Age*ProbeDuration, type="III")
  dataset.summary.3 = summary(dataset.anova.3, multivariate=FALSE)
  dataset.anovasheet.3 = as.matrix(dataset.summary.3$univariate.tests)
  
  dataset.4 = complete(dataset.mi, 4)
  data.lm.4 = lm(as.matrix(dataset.4[4:length(colnames(dataset.4))]) ~ 1 + Site * Genotype * Gender, data=dataset.4)
  dataset.anova.4 = Anova(data.lm.4, idata=idata,idesign=~Age*ProbeDuration, type="III")
  dataset.summary.4 = summary(dataset.anova.4, multivariate=FALSE)
  dataset.anovasheet.4 = as.matrix(dataset.summary.4$univariate.tests)
  
  dataset.5 = complete(dataset.mi, 5)
  data.lm.5 = lm(as.matrix(dataset.5[4:length(colnames(dataset.5))]) ~ 1 + Site * Genotype * Gender, data=dataset.5)
  dataset.anova.5 = Anova(data.lm.5, idata=idata,idesign=~Age*ProbeDuration, type="III")
  dataset.summary.5 = summary(dataset.anova.5, multivariate=FALSE)
  dataset.anovasheet.5 = as.matrix(dataset.summary.5$univariate.tests)
  
  dataset.list = list(data.lm.1,data.lm.2,data.lm.3,data.lm.4,data.lm.5)

  mi.anovafile = data.frame(matrix(nrow=32,ncol=7))
  colnames(mi.anovafile) = c('Test', "SS-III", 'num df', 'SSE-III', 'den df', 'F', 'p')
  mi.anovafile$Test = rownames(dataset.anovasheet.1)
  mi.anovafile$`num df` = as.vector(dataset.anovasheet.1[ ,2])
  mi.anovafile$`den df` = as.vector(dataset.anovasheet.1[ ,4])
  
  mi.ss.merge = data.frame(matrix(nrow=32, ncol=6))
  mi.ss.merge[ ,1] = as.vector(dataset.anovasheet.1[ ,1])
  mi.ss.merge[ ,2] = as.vector(dataset.anovasheet.2[ ,1])
  mi.ss.merge[ ,3] = as.vector(dataset.anovasheet.3[ ,1])
  mi.ss.merge[ ,4] = as.vector(dataset.anovasheet.4[ ,1])
  mi.ss.merge[ ,5] = as.vector(dataset.anovasheet.5[ ,1])
  for(a in 1:nrow(mi.ss.merge)){
    raw.values = c(mi.ss.merge[a,1], mi.ss.merge[a,2], mi.ss.merge[a,3], mi.ss.merge[a,4], mi.ss.merge[a,5])
    mi.ss.merge[a,6] = mean(raw.values)
  mi.anovafile$`SS-III` = as.vector(mi.ss.merge[ ,6])
  }
  mi.sse.merge = data.frame(matrix(nrow=32, ncol=6))
  mi.sse.merge[ ,1] = as.vector(dataset.anovasheet.1[ ,3])
  mi.sse.merge[ ,2] = as.vector(dataset.anovasheet.2[ ,3])
  mi.sse.merge[ ,3] = as.vector(dataset.anovasheet.3[ ,3])
  mi.sse.merge[ ,4] = as.vector(dataset.anovasheet.4[ ,3])
  mi.sse.merge[ ,5] = as.vector(dataset.anovasheet.5[ ,3])
  for(a in 1:nrow(mi.sse.merge)){
    raw.values = c(mi.sse.merge[a,1], mi.sse.merge[a,2], mi.sse.merge[a,3], mi.sse.merge[a,4], mi.sse.merge[a,5])
    mi.sse.merge[a,6] = mean(raw.values)
  }
  mi.anovafile$`SSE-III` = as.vector(mi.sse.merge[ ,6])
  
  for(a in 1:nrow(mi.anovafile)){
    temp.MST = mi.anovafile[a,2] / mi.anovafile[a,3]
    temp.MSE = mi.anovafile[a,4] / mi.anovafile[a,5]
    mi.anovafile[a,6] = temp.MST / temp.MSE
  }
  for(a in 1:8){
    mi.anovafile[a,7] = pf(mi.anovafile[a,6], mi.anovafile[a,3], mi.anovafile[a,5], lower.tail = FALSE)
  }
  for(a in 9:nrow(mi.anovafile)){
    #temp.epsi = 
    mi.anovafile[a,7] = pf(mi.anovafile[a,6], mi.anovafile[a,3] , mi.anovafile[a,5], lower.tail = FALSE)
  }
  
  return(mi.anovafile)
}



Data.GenerateMI.Function2 <- function(dataset){
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Sex')
  dataset.melt = melt(dataset, id.vars = c('AnimalID','Site','Strain','Genotype','Sex'))
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

Data.GenerateMM.Function <- function(dataset){
  colnames(dataset)[c(1:5)] = c('AnimalID','Site','Strain','Genotype','Gender')
  dataset.melt = melt(dataset, id.vars = c('AnimalID','Site','Strain','Genotype','Gender'))
  measure.list = as.data.frame(stringr::str_split_fixed(dataset.melt$variable, "_", 3))
  dataset.melt$Age = measure.list$V2
  dataset.melt$Age = as.factor(dataset.melt$Age)
  dataset.melt$ProbeDur = measure.list$V3
  dataset.melt$ProbeDur = as.factor(dataset.melt$ProbeDur)
  dataset.fixed = dataset.melt[ ,c(1,2,4,5,8,9,7)]
  colnames(dataset.fixed) = c('AnimalID','Site','Genotype','Gender','Age','ProbeDuration','Measure')
  #data.lmer = lmer(Measure ~ Site * Genotype * Gender * Age * ProbeDuration + (1 |AnimalID) + (1 |Age:AnimalID) + (1 | ProbeDuration:AnimalID), data=dataset, REML=FALSE)
  data.lme = lme(Measure ~ Site * Genotype * Gender * Age * ProbeDuration, random=list(AnimalID=pdBlocked(list(~1, pdIdent(~Age-1), pdIdent(~ProbeDuration-1)))), method="ML",data=dataset.fixed)
  data.anova = anova(data.lme)
  return(data.anova)
  }

## Pass ANOVA function through lists ##
data.3x.list.lmanova = lapply(data.3x.list.melt, Data.GenerateLM.Function, idata=data.3x.idata)
data.5x.list.lmanova = lapply(data.5x.list.melt, Data.GenerateLM.Function, idata=data.5x.idata)
data.app.list.lmanova = lapply(data.app.list.melt, Data.GenerateLM.Function, idata=data.app.idata)

#data.3x.list.mimputate = lapply(data.3x.list.melt, Data.GenerateMI.Function, idata=data.3x.idata)
#data.5x.list.mimputate = lapply(data.5x.list.melt, Data.GenerateMI.Function, idata=data.5x.idata)
#data.app.list.mimputate = lapply(data.app.list.melt, Data.GenerateMI.Function, idata=data.app.idata)

#data.3x.list.mixedmodelML = lapply(data.3x.list.melt, Data.GenerateMM.Function)
#data.5x.list.mixedmodelML = lapply(data.5x.list.probe, Data.GenerateMM.Function)
#data.app.list.mixedmodelML = lapply(data.app.list.probe, Data.GenerateMM.Function)

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


SME.Function = function(dataset, idata, measures, within.levels){
  SME.list = list()
  repeated.measure = TRUE
  both.repeated.present = FALSE
  dataset = dataset[complete.cases(dataset), ]
  idesign = "~Age*ProbeDuration"
  if(((is.element('Age',measures) == TRUE) & (is.element('ProbeDuration',measures) == FALSE) & (is.element('ProbeDuration',within.levels) == FALSE) & (is.element('Age',within.levels) == FALSE))){
    filter.list = list(6:9,10:13,14:17)
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=8))
    new.data.col[ ,c(1:5)] = dataset[ ,c(1:5)]
    new.data.col[ ,6] = apply(dataset[unlist(filter.list[1])], 1, mean, na.rm=TRUE)
    new.data.col[ ,7] = apply(dataset[unlist(filter.list[2])], 1, mean, na.rm=TRUE)
    new.data.col[ ,8] = apply(dataset[unlist(filter.list[3])], 1, mean, na.rm=TRUE)
    #if(is.element('APPPS1',dataset$Strain) == TRUE){
      #new.data.col[ ,9] = apply(dataset[ ,18:21], 1, mean, na.rm=TRUE)
    #}
    dataset = new.data.col
    colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','4Months','7Months','10Months')
    idesign = NA
    idata.edit = factor(c('4Months','7Months','10Months'))
    idata.edit = as.data.frame(idata.edit)
    #if(is.element('APPPS1',dataset$Strain) == TRUE){
      #colnames(dataset)[9] = "13Months"
      #idata.edit = factor(c('4Months','7Months','10Months', "13Months"))
      #idata.edit = as.data.frame(idata.edit)
    #}
    #colnames(idata.edit) = "Age"
    #idesign = as.formula(idesign)
  }
  if(((is.element('Age',measures) == FALSE) & (is.element('ProbeDuration',measures) == TRUE) & (is.element('ProbeDuration',within.levels) == FALSE) & (is.element('Age',within.levels) == FALSE))){
    filter.list = list(c(6,10,14), c(7,11,15), c(8,12,16), c(9,13,17))
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=9))
    new.data.col[ ,c(1:5)] = dataset[ ,c(1:5)]
    new.data.col[ ,6] = apply(dataset[unlist(filter.list[1])], 1, mean, na.rm=TRUE)
    new.data.col[ ,7] = apply(dataset[unlist(filter.list[2])], 1, mean, na.rm=TRUE)
    new.data.col[ ,8] = apply(dataset[unlist(filter.list[3])], 1, mean, na.rm=TRUE)
    new.data.col[ ,9] = apply(dataset[unlist(filter.list[4])], 1, mean, na.rm=TRUE)
    dataset = new.data.col
    colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','600ms','800ms','1000ms','1500ms')
    idata.edit = factor(c('0600ms','0800ms','1000ms','1500ms'))
    idata.edit = as.data.frame(idata.edit)
    #colnames(idata.edit) = "ProbeDuration"
    #idesign = NA
    #idesign = as.formula(idesign)
  }
  if(((is.element('Age',measures) == FALSE) & (is.element('ProbeDuration',measures) == FALSE) & (is.element('ProbeDuration',within.levels) == FALSE) & (is.element('Age',within.levels) == TRUE))){
    filter.list = list(6:9,10:13,14:17)
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=8))
    new.data.col[ ,c(1:5)] = dataset[ ,c(1:5)]
    new.data.col[ ,6] = apply(dataset[unlist(filter.list[1])], 1, mean, na.rm=TRUE)
    new.data.col[ ,7] = apply(dataset[unlist(filter.list[2])], 1, mean, na.rm=TRUE)
    new.data.col[ ,8] = apply(dataset[unlist(filter.list[3])], 1, mean, na.rm=TRUE)
    #if(is.element('APPPS1',dataset$Strain) == TRUE){
      #new.data.col[ ,9] = apply(dataset[ ,18:21], 1, mean, na.rm=TRUE)
    #}
    dataset = new.data.col
    idesign = NA
    colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','04','07','10')
    #if(is.element('APPPS1',dataset$Strain) == TRUE){
      #colnames(dataset)[9] = "13Months"
    #}
    repeated.measure = FALSE
  }
  if(((is.element('Age',measures) == FALSE) & (is.element('ProbeDuration',measures) == FALSE) & (is.element('ProbeDuration',within.levels) == TRUE) & (is.element('Age',within.levels) == FALSE))){
    filter.list = list(c(6,10,14), c(7,11,15), c(8,12,16), c(9,13,17))
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=9))
    new.data.col[ ,c(1:5)] = dataset[ ,c(1:5)]
    new.data.col[ ,6] = apply(dataset[ ,unlist(filter.list[1])], 1, mean, na.rm=TRUE)
    new.data.col[ ,7] = apply(dataset[unlist(filter.list[2])], 1, mean, na.rm=TRUE)
    new.data.col[ ,8] = apply(dataset[unlist(filter.list[3])], 1, mean, na.rm=TRUE)
    new.data.col[ ,9] = apply(dataset[unlist(filter.list[4])], 1, mean, na.rm=TRUE)
    dataset = new.data.col
    colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','600ms','800ms','1000ms','1500ms')
    idesign = NA
    repeated.measure = FALSE
  }
    
  if(((is.element('Age',measures) == FALSE) & (is.element('ProbeDuration',measures) == FALSE) & (is.element('Age',within.levels) == FALSE) & (is.element('ProbeDuration',within.levels) == FALSE))){
    new.data.col = as.data.frame(matrix(nrow = nrow(dataset),ncol=6))
    new.data.col[ ,c(1:5)] = dataset[ ,c(1:5)]
    new.data.col[ ,6] = apply(dataset[6:17], 1, mean, na.rm = TRUE)
    dataset = new.data.col
    colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','Measure')
    idesign = NA
    repeated.measure = FALSE
  }
  if(((is.element('Age',measures) == TRUE) & (is.element('ProbeDuration',measures) == TRUE) & (is.element('Age',within.levels) == FALSE) & (is.element('ProbeDuration',within.levels) == FALSE))){
    idesign = "~Age*ProbeDuration"
    idesign = as.formula(idesign)
    idata.edit = idata
    both.repeated.present = TRUE
  }
  if(isTRUE(within.levels == "Gender")){
    filter.list = list('M',"F")
    filter.col = 5
  }else if(isTRUE(within.levels == "Genotype")){
    filter.list = list('w','t')
    filter.col = 4
  }else if(isTRUE(within.levels == "Site")){
    filter.list = list('UOG', 'UWO')
    filter.col = 2
  }else if(isTRUE(within.levels == "Age")){
    if(is.element('ProbeDuration',measures) == TRUE){
      filter.list = list(6:9,10:13,14:17)
      #if(is.element('APPPS1',dataset$Strain) == TRUE){
        #filter.list = list(6:9,10:13,14:17,18:21)
      #}
      idata.edit = factor(c('0600ms','0800ms','1000ms','1500ms'))
      idata.edit = as.data.frame(idata.edit)
      #colnames(idata.edit) = "ProbeDuration"
      idesign = "~ProbeDuration"
      idesign = as.formula(idesign)
    }else{
      filter.list = list(6,7,8)
      #if(is.element('APPPS1',dataset$Strain) == TRUE){
        #filter.list = list(6,7,8,9)
      #}
      idata.edit = NA
      idesign = NA
    }
    filter.col = NA
  }else if(isTRUE(within.levels == "ProbeDuration")){
    if(is.element('Age',measures) == TRUE){
      filter.list = list(c(6,10,14), c(7,11,15), c(8,12,16), c(9,13,17))
      idata.edit = factor(c('4Months','7Months','10Months'))
      #if(is.element('APPPS1',dataset$Strain) == TRUE){
        #filter.list = list(c(6,10,14,18), c(7,11,15,19), c(8,12,16,20), c(9,13,17,21))
        #idata.edit = factor(c('4Months','7Months','10Months','13Months'))
      #}
      idata.edit = as.data.frame(idata.edit)
      #colnames(idata.edit) = "Age"
      idesign = "~Age"
      idesign = as.formula(idesign)
    }else{
      filter.list = list(6,7,8,9)
      idata.edit = NA
      idesign= NA
    }
    filter.col = NA
  }
  if(is.element('Age',measures) == TRUE){
    measures = measures[measures != "Age"]
  }
  if(is.element('ProbeDuration',measures) == TRUE){
    measures = measures[measures != "ProbeDuration"]
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
      dataset.anova = dataset[ ,c(1:5, current.filter)]
      dataset.anova = dataset.anova[complete.cases(dataset.anova), ]
      data.matrix = as.matrix(dataset.anova[ ,6:ncol(dataset.anova)])
      data.lm = lm(measure.string, data = dataset.anova)
      data.anova = Anova(data.lm, idata=idata.edit,idesign= ~idata.edit, type="III")
      data.summary = summary(data.anova, multivariate=FALSE)
      #data.summary = data.summary$univariate.tests
      SME.list[[a]] = data.summary
    }
  }else if((is.na(filter.col) == FALSE) & (repeated.measure == TRUE) & (both.repeated.present == FALSE) ){
    for(a in 1:length(filter.list)){
      dataset.anova = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ]
      dataset.anova = dataset.anova[complete.cases(dataset.anova), ]
      data.matrix = as.matrix(dataset.anova[, 6:ncol(dataset.anova)])
      data.lm = lm(measure.string, data = dataset.anova)
      data.anova = Anova(data.lm, idata=as.data.frame(idata.edit),idesign=~idata.edit, type="III")
      data.summary = summary(data.anova, multivariate=FALSE)
      #data.summary = data.summary$univariate.tests
      SME.list[[a]] = data.summary
    }
  }else if((is.na(filter.col) == FALSE) & (repeated.measure == TRUE) & (both.repeated.present == TRUE)){
    for(a in 1:length(filter.list)){
      dataset.anova = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ]
      dataset.anova = dataset.anova[complete.cases(dataset.anova), ]
      data.matrix = as.matrix(dataset.anova[, 6:ncol(dataset.anova)])
      data.lm = lm(measure.string, data = dataset.anova)
      data.anova = Anova(data.lm, idata=idata.edit,idesign=idesign, type="III")
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
      data.matrix = as.matrix(dataset[which(dataset[filter.col] == unlist(filter.list[a])), 6:ncol(dataset)])
      data.lm = lm(measure.string, data = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ])
      data.anova = Anova(data.lm, type="III")
      SME.list[[a]] = data.anova
    }
  }
  return(SME.list)
}
  
test2 = SME.Function(data.app.list.melt$data.app.probe.omission, data.app.idata, c('Age'), 'Site')

SME.Function.2 = function(dataset, between.measure,within.measure, within.levels){
  dataset.complete = dataset[complete.cases(dataset), ]
  dataset.2 = melt(dataset.complete, id.vars=c('AnimalID','Site','Strain','Genotype','Gender'))
  dataset.within = as.data.frame(stringr:: str_split_fixed(dataset.2$variable,"_",2))
  dataset = cbind(dataset.2,dataset.within)
  dataset = dataset[ ,c(1:5,8:9,7)]
  colnames(dataset) = c('AnimalID','Site','Strain','Genotype','Gender','Age','ProbeDuration','Measure')
  
  SME.list = list()
  
  if(isTRUE(within.levels == "Gender")){
    filter.list = list('M',"F")
    filter.col = 5
  }else if(isTRUE(within.levels == "Genotype")){
    filter.list = list('w','t')
    filter.col = 4
  }else if(isTRUE(within.levels == "Site")){
    filter.list = list('UOG', 'UWO')
    filter.col = 2
  }else if(isTRUE(within.levels == "Age")){
    filter.list = list('04', '07', '10')
    filter.col = 6
  }else if(isTRUE(within.levels == "ProbeDuration")){
    filter.list = list('0600ms', '0800ms', '1000ms', '1500ms')
    filter.col = 7
  }
  ez.anova.text.1 = "ezANOVA("
  for(a in 1:length(filter.list)){
    filter.data = dataset[which(dataset[filter.col] == unlist(filter.list[a])), ]
    if(is.null(between.measure) == TRUE){
      anova.data = ezANOVA(data=filter.data, wid= AnimalID, dv= Measure, within = parse(text=within.measure) )
    }else if(is.null(within.measure) == TRUE){
      anova.data = ezANOVA(data=filter.data, wid= AnimalID, dv= Measure, between = eval(between.measure) )
    }else{
      anova.data = ezANOVA(data=filter.data, wid= AnimalID, dv= Measure, between = eval(between.measure), within = eval(within.measure) )
    }
    SME.list[[a]] = anova.data
  }
  return(SME.list)
}

SME.Function(data.app.list.melt$data.app.probe.persev , idata= data.app.idata, measure = c('Site','Gender','ProbeDuration'), within.levels = 'Age')

