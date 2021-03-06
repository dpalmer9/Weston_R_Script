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

## Function Version ##
SME.Calculation.Function = function(newSST, newDF1, oldSSEWTH,oldSSEBTW=NULL,oldDFWTH,oldDFBTW=NULL, within.type,factor.types=NULL){
  within.type.btw = c('Between', "Btw", "between", '1',1)
  within.type.wth = c('Within', 'Wth', "within", '2',2)
  if((is.element(within.type, within.type.btw) == TRUE) & (isTRUE(factor.types != 1))){
    SME.error = oldSSEWTH
    SME.df2 = oldDFWTH
    SME.df1 = newDF1
    SME.SST = newSST
    SME.F = (SME.SST / SME.df1) / (SME.error / SME.df2)
    SME.Fsig = pf(SME.F, SME.df1, SME.df2, lower.tail = FALSE)
  }else if((is.element(within.type, within.type.btw) == TRUE) & (isTRUE(factor.types == 1))){
    SME.error = oldSSEBTW
    SME.df2 = oldDFBTW
    SME.df1 = newDF1
    SME.SST = newSST
    SME.F = (SME.SST / SME.df1) / (SME.error / SME.df2)
    SME.Fsig = pf(SME.F, SME.df1, SME.df2, lower.tail = FALSE)
  }else if(is.element(within.type, within.type.wth) == TRUE){
    SME.error = (oldSSEWTH + oldSSEBTW) / (oldDFWTH + oldDFBTW)
    SME.df2 = ((oldSSEWTH + oldSSEBTW) ^ 2) / (((oldSSEWTH ^2) / oldDFWTH) + ((oldSSEBTW ^2) / oldDFBTW))
    SME.df1 = newDF1
    SME.SST = newSST
    SME.F = (SME.SST / SME.df1) / SME.error
    SME.Fsig = pf(SME.F, SME.df1, SME.df2, lower.tail = FALSE)
  }
  final.data = as.data.frame(matrix(nrow=1,ncol=6))
  final.data[1, ] = c(SME.SST,SME.error,SME.df1,SME.df2,SME.F,SME.Fsig)
  colnames(final.data) = c('SST','SSE','df1','df2','F','sig(p)')
  return(final.data)
}
# Between #
SME.Calculation.Function(newSST = 0.1402, newDF1 = 1, oldSSEWTH = 56741, oldDFWTH = 396, oldSSEBTW = 65.038 ,oldDFBTW = 50 ,within.type = "Between", factor.types=1)
# Within #
SME.Calculation.Function(newSST = 26.82, newDF1 = 1, oldSSEWTH = 618.52, oldDFWTH = 495, oldSSEBTW = 506.6 ,oldDFBTW = 55 , within.type = "Within")


## Boniferroni Correction ##
Boniferonni.Func <- function(pvalue){
  alpha1 = 0.05
  alpha2 = 0.01
  alpha3 = 0.001
}

p.adjust(p= 0.000042226 ,method = "bonferroni", n= 2)
