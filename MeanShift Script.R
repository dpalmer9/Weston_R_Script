## Mean Shift Algorithm ##

## Factors ##
set.seed(100)
ms.iterations = 1000
ms.bandwidth = c(0.5,0.5)


## Read Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

raw.data.probe = raw.data.probe[ ,c(3,4,5,6,7,8,9,16,17,18,19,122,70)]
colnames(raw.data.probe) = c('Animal.ID','TestSite','MouseStrain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature','Perseverative','CorrectLatency','RewardLatency')

## Omit Data ##
raw.data.probe = na.omit(raw.data.probe)

## Scale Data ##
raw.data.probe[ ,c(8:13)] = scale(raw.data.probe[ ,c(8:13)])

train.data = as.matrix(raw.data.probe[ ,c(8:13)])


## Get Bandwith ##

bandwith.check = quantile(dist(t(train.data)),seq(0.05,0.40,by=0.05))

## Mean Shift ##
mean.shift.data = meanShiftR::meanShift(train.data,train.data
                                        ,alpha=0,iterations = ms.iterations)
