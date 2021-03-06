## Library ##
library(kml3d)
library(reshape2)

## Functions ##
Data.Formatting.Function.5CSRTT = function(dataset){
  dataset$StimulusLength = as.factor(dataset$StimulusLength)
  dataset$StimulusLength = relevel(dataset$StimulusLength,'0.8')
  dataset$StimulusLength = relevel(dataset$StimulusLength,'1')
  dataset$StimulusLength = relevel(dataset$StimulusLength,'1.5')
  melt.data = melt(dataset,id.vars=c(1:7))
  cast.data = dcast(melt.data,Animal.ID + TestSite + MouseStrain + Genotype + Sex ~ variable + Age + StimulusLength, fun.aggregate = mean, na.rm=TRUE)
  return(cast.data)
}


## Read Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

raw.data.probe = raw.data.probe[ ,c(3,4,5,6,7,8,9,16,17,18,19,122,70)]
colnames(raw.data.probe) = c('Animal.ID','TestSite','MouseStrain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature','Perseverative','CorrectLatency','RewardLatency')

raw.data.probe[, c(8:13)] = scale(raw.data.probe[, c(8:13)])

formatted.data.probe = Data.Formatting.Function.5CSRTT(raw.data.probe)

formatted.data.probe = na.omit(formatted.data.probe) #Omit Data


## Set Parameters ##
kml.group.no = 3

## Rum K-Mean Algorithm ##

kma.3d.analysis = cld3d(formatted.data.probe, timeInData = list(Accuracy=c(6:17)
                                                                ,Omission=c(18:29)
                                                                ,Premature=c(30:41)
                                                                ,Perseverative=c(42:53)
                                                                ,CorrectLatency=c(54:65)
                                                                ,RewardLatency=c(66:77)))

kml3d(kma.3d.analysis)

combined.data = formatted.data.probe

combined.data$Cluster = getClusters(kma.3d.analysis,kml.group.no)


## Visualize
graphing.data.mean = aggregate(combined.data[ ,6:77],by=list(combined.data$Cluster), FUN=mean, na.rm=TRUE)

graphing.melt.data = melt(graphing.data.mean)

new.group.data = as.data.frame(stringr::str_split_fixed(graphing.melt.data$variable,'_',2))

day.vector = as.vector(new.group.data[ ,2])
day.unique = as.vector(unique(day.vector))

for(a in 1:length(day.unique)){
  day.vector[day.vector==day.unique[a]] = a
}

graphing.melt.data[ ,2] = new.group.data[ ,1]
graphing.melt.data$Time = day.vector

graphing.data.cast = dcast(graphing.melt.data, Group.1 + Time ~ variable, fun.aggregate = mean, na.rm=TRUE)
graphing.data.cast = graphing.data.cast[order(graphing.data.cast$Group.1,as.numeric(graphing.data.cast$Time)), ]

graphing.data.cast.4m = graphing.data.cast[c(1:4,13:16,25:28), ]
graphing.data.cast.4m$Time = dplyr::recode_factor(graphing.data.cast.4m$Time, '4'='0.6', '3'='0.8', '2'='1', '1'='1.5')
graphing.data.cast.4m$Time = factor(graphing.data.cast.4m$Time,levels=graphing.data.cast.4m$Time[c(4,3,2,1)])
graphing.data.cast.7m = graphing.data.cast[c(5:8,17:20,29:32), ]
graphing.data.cast.7m$Time = dplyr::recode_factor(graphing.data.cast.7m$Time, '8'='0.6', '7'='0.8', '6'='1', '5'='1.5')
graphing.data.cast.7m$Time = factor(graphing.data.cast.7m$Time,levels=graphing.data.cast.7m$Time[c(4,3,2,1)])
graphing.data.cast.10m = graphing.data.cast[c(9:12,21:24,33:36), ]
graphing.data.cast.10m$Time = dplyr::recode_factor(graphing.data.cast.10m$Time, '12'='0.6', '11'='0.8', '10'='1', '9'='1.5')
graphing.data.cast.10m$Time = factor(graphing.data.cast.10m$Time,levels=graphing.data.cast.10m$Time[c(4,3,2,1)])

ggplot(data=graphing.data.cast.4m,aes(x=graphing.data.cast.4m$Time,y=graphing.data.cast.4m$Accuracy,group=as.factor(graphing.data.cast.4m$Group.1))) + 
  geom_line(aes(color=as.factor(graphing.data.cast.4m$Group.1))) + 
  geom_point(aes(color=as.factor(graphing.data.cast.4m$Group.1))) +
  scale_x_discrete(limits = rev(levels(graphing.data.cast.4m$Time)))

ggplot(data=graphing.data.cast.7m,aes(x=as.factor(graphing.data.cast.7m$Time),y=graphing.data.cast.7m$Accuracy,group=as.factor(graphing.data.cast.7m$Group.1))) + 
  geom_line(aes(color=as.factor(graphing.data.cast.7m$Group.1))) + 
  geom_point(aes(color=as.factor(graphing.data.cast.7m$Group.1)))

ggplot(data=graphing.data.cast.10m,aes(x=as.factor(graphing.data.cast.10m$Time),y=graphing.data.cast.10m$Accuracy,group=as.factor(graphing.data.cast.10m$Group.1))) + 
  geom_line(aes(color=as.factor(graphing.data.cast.10m$Group.1))) + 
  geom_point(aes(color=as.factor(graphing.data.cast.10m$Group.1)))
