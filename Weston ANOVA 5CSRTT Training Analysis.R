## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
library(plyr)


## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceAcq.csv')
colnames(raw.data) = c('AnimalID','Site','Strain','Genotype','Gender','Duration')

## Separate By Strain ##
data.3x = raw.data[which(raw.data$Strain==" 3XTG"), ]
data.3x$Strain = NULL
data.5x = raw.data[which(raw.data$Strain==" 5FAD"), ]
data.5x$Strain = NULL
data.app = raw.data[which(raw.data$Strain==" APPPS1"), ]
data.app$Strain = NULL

## Separate by Probe Duration & Generate Frequency ##
Data.ProbeSeparation.Function <- function(dataset){
  datalist = list()
  datalist$data.four = dataset[which(dataset$Duration == "4"), ]
  datalist$data.four$Duration = NULL
  datalist$data.four = plyr::count(datalist$data.four, c('AnimalID','Site','Genotype','Gender'))
  datalist$data.two = dataset[which(dataset$Duration == "2"), ]
  datalist$data.two$Duration = NULL
  datalist$data.two = plyr::count(datalist$data.two, c('AnimalID','Site','Genotype','Gender'))
  return(datalist)
}
data.3x.list = Data.ProbeSeparation.Function(data.3x)
data.5x.list = Data.ProbeSeparation.Function(data.5x)
data.app.list = Data.ProbeSeparation.Function(data.app)

## Generate ANOVA ##
Data.TrainingANOVA.Function <- function(datafile){
  temp.lm = lm(freq ~ Site * Genotype * Gender, data=datafile)
  anova.result = Anova(temp.lm, type="III")
  return(anova.result)
}

data.3x.anova = lapply(data.3x.list, Data.TrainingANOVA.Function)
data.5x.anova = lapply(data.5x.list, Data.TrainingANOVA.Function)
data.app.anova = lapply(data.app.list, Data.TrainingANOVA.Function)
