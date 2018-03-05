PAL_Data <- read.csv("~/Documents/PAL_Data.csv")# You will have to change this line for the appropriate file.#
PAL_Data = na.omit(PD_Data)
names(PAL_Data) # This will be useful to create the data matrix

k=PD_Data
m=as.data.frame(k)
#Convert Matrix to cld
library("kml3d")
#Create Data Matrix you will need to edit the c(x,y,z) in order for it to match the appropriate column number which you can find using the names() function#

cld3dPregTemp <- cld3d(m,timeInData = list(Errors=c(8,15,22,29,36,43,50,57,64)
                                           ,Accuracy=c(9,16,23,30,37,44,51,58,65)
                                           ,CompTime=c(6,13,20,27,34,41,48,55,62)
                                           ,TrialsCompleted=c(7,14,21,28,35,42,49,56,63)
                                           ,Touch=c(10,17,24,31,38,45,52,59,66)
                                           ,Reward=c(11,18,25,32,39,46,53,60,67)))

kml3d(cld3dPregTemp)
plotMeans3d(cld3dPregTemp,3)

m$clusters <- getClusters(cld3dPregTemp, 3)

#This will will write a new file which is identical to the input file, but with a new column at the very end with the cluster that mouse belongs to. You will then have to do the z-scores of each variable used in the analysis and then re-organize the cluster names. It will call them A,B,C but in no order.# 
write.csv(m, file="PAL_ClusterResults.csv")
