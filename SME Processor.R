## Options ##
options(scipen=100)
options(contrasts = c('contr.sum','contr.poly'))

## Library ##
library(tidyverse)

## Get Input
SST.value = readline(prompt = "Please input SST from lower analysis:")
dfn.value = readline(prompt = "Please input first degree of freedom:")
SSE.value = readline(prompt = "Please input SSE from higher analysis:")
dfd.value = readline(prompt = "Please input second degree of freedom:")

SST.value = as.numeric(SST.value)
dfn.value = as.numeric(dfn.value)
SSE.value = as.numeric(SSE.value)
dfd.value = as.numeric(dfd.value)
## Calculate ##
new.MSE = (SSE.within.value + SSE.between.value) / (df.within.value + df.between.value)
new.df = ((SSE.within.value + SSE.between.value)^2) / (((SSE.within.value)^2 / df.within.value) + ((SSE.between.value)^2 / df.between.value))
new.F = (SST.value / dfn.value) / new.MSE
F.sig = pf(new.F, dfn.value, new.df)
