## Library ##
library(tidyverse)
library(reshape2)
library(car)
library(lme4)
library(nlm3)

## Acquire data ##
raw.data = read.csv('C:\\Users\\dpalmer\\Documents\\WestonANOVAProcedure\\PD Reversal.csv')

## Re Calculate Mean Value for Error ##


## Remove Individual Latency ##
raw.data.clean = raw.data[ ,c(2:8,10,13:16,)]


## Separate By Strain ##
data.3x = raw.data[which(raw.data$Mouse.strain==" 3XTG"), ]
data.5x = raw.data[which(raw.data$Mouse.strain==" 5FAD"), ]
data.app = raw.data[which(raw.data$Mouse.strain==" APPPS1"), ]

