## Library ##
library(tidyverse)
library(lme4)

## Parameters ##
set.seed(100)

## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\FiveChoiceProbe.csv')

data.3x = raw.data[which(raw.data$Mouse.strain=="3XTG"), ]
data.5x = raw.data[which(raw.data$Mouse.strain=="5FAD"), ]
data.app = raw.data[which(raw.data$Mouse.strain=="APPPS1"), ]
## Format Data ##


## Generate Lm ##
lmer.3x.acc = lmer(Threshold...Accuracy..~ Genotype + Gender + Age..months. + Genotype*Gender + Genotype*Age..months. + Gender*Age..months. + Genotype*Gender*Age..months. + (1|Animal.ID), data=data.3x,REML=TRUE)
lmer.5x.acc = lmer(Threshold...Accuracy..~ Genotype + Gender + Age..months. + Genotype*Gender + Genotype*Age..months. + Gender*Age..months. + Genotype*Gender*Age..months. + (1|Animal.ID), data=data.5x,REML=TRUE)
lmer.app.acc = lmer(Threshold...Accuracy..~ Genotype + Gender + Age..months. + Genotype*Gender + Genotype*Age..months. + Gender*Age..months. + Genotype*Gender*Age..months. + (1|Animal.ID),data=data.app,REML=TRUE)