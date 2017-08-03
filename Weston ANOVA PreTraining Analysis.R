## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
#library(nlm3)

## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoicePreTrain.csv')
colnames(raw.data) = c('AnimalID','Site','Strain','Genotype','Gender','Age','Protocol')
raw.data$Age = NULL

## Switch Protocol Names to More Meaningful Label ##
Data.ProtocolNameSwap.Function <- function(dataset){
  dataset$Protocol = as.character(dataset$Protocol)
  dataset$Protocol[dataset$Protocol == "5CSRTT_Habituation_1_"] = "01_Habituation1"
  dataset$Protocol[dataset$Protocol == "5CSRTT_Habituation_2_"] = "02_Habituation2"
  dataset$Protocol[dataset$Protocol == "5CSRTT_Initial_Touch_Training_v3"] = "03_InitialTouch"
  dataset$Protocol[dataset$Protocol == " 5CSRTT_Initial_Touch_Training_v3"] = "03_InitialTouch"
  dataset$Protocol[dataset$Protocol == "5CSRTT_Must_Touch_Training_v2"] = "04_MustTouch"
  dataset$Protocol[dataset$Protocol == " 5CSRTT_Must_Touch_Training_v2"] = "04_MustTouch"
  dataset$Protocol[dataset$Protocol == "5CSRTT_Must_Initiate_Training_v1"] = "05_MustInitiate"
  dataset$Protocol[dataset$Protocol == "5CSRTT_Punish_Incorrect_Training_v3"] = "06_PunishIncorrect"
  dataset$Protocol[dataset$Protocol == " 5CSRTT_Punish_Incorrect_Training_v3"] = "06_PunishIncorrect"
  dataset$Protocol = as.factor(dataset$Protocol)
  return(dataset)
}

raw.data = Data.ProtocolNameSwap.Function(raw.data)
## Separate By Strain ##
data.3x = raw.data[which(raw.data$Strain=="3xTG"), ]
data.3x$Strain = NULL
data.5x = raw.data[which(raw.data$Strain=="5xFAD"), ]
data.5x$Strain = NULL
data.app = raw.data[which(raw.data$Strain=="APPPS1"), ]
data.app$Strain = NULL

## Generate Freq ##
data.3x.freq = plyr::count(data.3x, c('AnimalID','Site','Genotype','Gender','Protocol'))
data.5x.freq = plyr::count(data.5x, c('AnimalID','Site','Genotype','Gender','Protocol'))
data.app.freq = plyr::count(data.app, c('AnimalID','Site','Genotype','Gender','Protocol'))

## Separate Sheets Function ##
Data.FreqSplit.Function <- function(dataset){
  datalist = list()
  datalist$data.initial = dataset[which(dataset$Protocol == "03_InitialTouch"), ]
  datalist$data.musttouch = dataset[which(dataset$Protocol == "04_MustTouch"), ]
  datalist$data.mustinitiate = dataset[which(dataset$Protocol == "05_MustInitiate"), ]
  datalist$data.punishincorrect = dataset[which(dataset$Protocol == "06_PunishIncorrect"), ]
  return(datalist)
}

data.3x.freq.list = Data.FreqSplit.Function(data.3x.freq)
data.5x.freq.list = Data.FreqSplit.Function(data.5x.freq)
data.app.freq.list = Data.FreqSplit.Function(data.app.freq)

## ANOVA Function ##
Data.ANOVA.Function <- function(datafile){
    datafile$Protocol = NULL
    temp.lm = lm(freq ~ Site * Genotype * Gender, data=datafile)
    anova.result = Anova(temp.lm, type="III")
    return(anova.result)
}


data.3x.anova = lapply(data.3x.freq.list, Data.ANOVA.Function)
data.5x.anova = lapply(data.5x.freq.list, Data.ANOVA.Function)
#app Initial Touch has no variability, remove from list
data.app.freq.list$data.initial = NULL
data.app.anova = lapply(data.app.freq.list, Data.ANOVA.Function)
