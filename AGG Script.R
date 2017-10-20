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
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceDatabase2Second.csv')
#raw.data = read.csv('C:\\Users\\Danie\\Documents\\R\\Projects\\Weston_R_Script\\Weston PAL Pretrain.csv')
raw.data$Age..months. = as.factor(raw.data$Age..months.)

#melt.data = melt(raw.data, id.vars = c("Database","Animal.ID","Site","Mouse.strain","GenoAPPPS1ype","Gender","Age..months."))

agg.list = list(raw.data$Animal.ID,raw.data$Site,raw.data$Mouse.strain,raw.data$Genotype,raw.data$Gender,raw.data$Age..months.,raw.data$Schedule.name)
raw.data$Max.number.of.trials = NULL
raw.data$Max.schedule.time = NULL
raw.data$Schedule.run.date = NULL
agg.data.2 = aggregate(raw.data.2, by= agg.list, FUN=mean, na.rm = TRUE)
agg.data.1 = aggregate(Day ~ Animal.ID + Site + Mouse.strain + Genotype + Gender + Age..months. + Schedule.name, FUN=sum, na.rm=TRUE, data=raw.data)

agg.data.2[ ,8:14] = agg.data.2 [ ,1:7]
agg.data.2[ ,1:7] = NULL

write.csv(agg.data.1, "Weston 5CSRTT 2 Second Training AGG.csv")

raw.data$Max.number.of.trials