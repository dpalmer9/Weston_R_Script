## Library ##
library(kml3d)
library(reshape2)

## Functions ##
Data.Formatting.Function.5CSRTT = function(dataset){
  melt.data = melt(dataset,id.vars=c(1:7))
  cast.data = dcast(melt.data,Animal.ID + TestSite + MouseStrain + Genotype + Sex + Age ~ variable + StimulusLength, fun.aggregate = mean, na.rm=TRUE)
  return(cast.data)
}


## Read Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

raw.data.probe = raw.data.probe[ ,c(3,4,5,6,7,8,9,16,17,18,19,122,70)]
colnames(raw.data.probe) = c('Animal.ID','TestSite','MouseStrain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature','Perseverative','CorrectLatency','RewardLatency')

formatted.data.probe = Data.Formatting.Function.5CSRTT(raw.data.probe)

formatted.data.probe = na.omit(formatted.data.probe) #Omit Data

## Set Parameters ##
kml.group.no = 3

## Rum K-Mean Algorithm ##

kma.3d.analysis = cld3d(formatted.data.probe, timeInData = list(Accuracy=c(10,9,8,7)
                                                                ,Omission=c(14,13,12,11)
                                                                ,Premature=c(18,17,16,15)
                                                                ,Perseverative=c(22,21,20,19)
                                                                ,CorrectLatency=c(26,25,24,23)
                                                                ,RewardLatency=c(30,29,28,27)))
