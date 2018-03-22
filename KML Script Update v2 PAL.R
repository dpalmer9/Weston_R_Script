## Library ##
library(dplyr)
library(kml3d)
library(reshape2)
library(ggplot2)

## Functions ##
Data.Formatting.Function.PAL = function(dataset){
  melt.data = melt(dataset,id.vars=c(1:7))
  melt.data[ ,1] = as.character(melt.data[ ,1])
  melt.data[ ,6] = as.character(melt.data[ ,6])
  for(a in 1:nrow(melt.data)){
    melt.data[a,1] = paste(as.character(melt.data[a,1]),as.character(melt.data[a,6]),sep='.')
  }
  cast.data = dcast(melt.data,Animal.ID + MouseStrain + Genotype + Sex + Age ~ variable + Week, fun.aggregate = mean, na.rm=TRUE)
  return(cast.data)
}


## Read Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PAL\\Weston PAL Main Task Aggregated Oct 12 2017 Updated.csv')

raw.data.probe = raw.data.probe[ ,c(2:7,10,15:16,239,315)]
colnames(raw.data.probe) = c('Animal.ID','TestSite','MouseStrain','Genotype','Sex','Age','Week','Corrections','Accuracy','CorrectLatency','RewardLatency')

for(a in 8:11){
  raw.data.probe[,a] = scale(raw.data.probe[,a])
}


formatted.data.probe = Data.Formatting.Function.PAL(raw.data.probe)

formatted.data.probe = na.omit(formatted.data.probe) #Omit Data


## Set Parameters ##
kml.group.no = 3

## Rum K-Mean Algorithm ##

kma.3d.analysis = cld3d(formatted.data.probe, timeInData = list(Corrections=c(6:14)
                                                                ,Accuracy=c(15:23)
                                                                ,CorrectLatency=c(24:32)
                                                                ,RewardLatency=c(33:41)))

kml3d(kma.3d.analysis)

combined.data = formatted.data.probe

combined.data$Cluster = getClusters(kma.3d.analysis,3)

## Save Count ##
count.data = combined.data
count.data$Count = 1
count.data = aggregate(count.data$Count,by=list(count.data$Cluster,count.data$Genotype,count.data$Sex,count.data$Age),FUN=sum, na.rm=TRUE)
colnames(count.data) = c('Cluster','Strain','Sex','Age','Count')
count.data$Cluster = recode(count.data$Cluster,'A' = 'High','B' = 'Low', 'C' = 'Mid')
## Visualize
graphing.data.mean = aggregate(combined.data[ ,6:41],by=list(combined.data$Cluster), FUN=mean, na.rm=TRUE)

graphing.melt.data = melt(graphing.data.mean)

new.group.data = as.data.frame(stringr::str_split_fixed(graphing.melt.data$variable,'_',2))
day.vector = new.group.data[ ,2]

graphing.melt.data$variable = new.group.data[ ,1]
graphing.melt.data$Time = day.vector

graphing.data.cast = dcast(graphing.melt.data, Group.1 + Time ~ variable, fun.aggregate = mean, na.rm=TRUE)
graphing.data.cast = graphing.data.cast[order(graphing.data.cast$Group.1,as.numeric(graphing.data.cast$Time)), ]

ggplot(data=graphing.data.cast,aes(x=graphing.data.cast$Time,y=graphing.data.cast$Accuracy,group=as.factor(graphing.data.cast$Group.1))) + 
  geom_line(aes(color=as.factor(graphing.data.cast$Group.1))) + 
  geom_point(aes(color=as.factor(graphing.data.cast$Group.1)))
