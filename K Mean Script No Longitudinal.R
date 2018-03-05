## K Mean Alternate ##

## Library ##

## Factors ##
kml.group.no = 3

## Read Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

raw.data.probe = raw.data.probe[ ,c(3,4,5,6,7,8,9,16,17,18,19,122,70)]
colnames(raw.data.probe) = c('Animal.ID','TestSite','MouseStrain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature','Perseverative','CorrectLatency','RewardLatency')

## Omit Data ##
raw.data.probe = na.omit(raw.data.probe)

## Scale Data ##
raw.data.probe[ ,c(8:13)] = scale(raw.data.probe[ ,c(8:13)])

## K-Mean (No Longitudinal): Ignore Factors-Unsupervised ##
kmean.analysis.data = kmeans(raw.data.probe[ ,c(8:13)],kml.group.no)


## Merge Data from K-Mean
kmean.data = raw.data.probe
kmean.data$KGroups = kmean.analysis.data$cluster

## Aggregate by K-Groups ##
kmean.data.agg = aggregate(kmean.data[ ,c(8:13)], by=list(kmean.data$KGroups, kmean.data$StimulusLength),FUN=mean,na.rm=TRUE)


## Visualize ##
ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$Accuracy,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))

ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$Omissions,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))

ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$Premature,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))

ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$Perseverative,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))

ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$CorrectLatency,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))

ggplot(data=kmean.data.agg,aes(x=as.factor(kmean.data.agg$Group.2),y=kmean.data.agg$RewardLatency,group=as.factor(kmean.data.agg$Group.1))) + 
  geom_line(aes(color=as.factor(kmean.data.agg$Group.1))) + 
  geom_point(aes(color=as.factor(kmean.data.agg$Group.1)))
