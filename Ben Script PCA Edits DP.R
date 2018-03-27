## Library ##
library(ggplot2)
library(ggfortify)
library(factoextra)

#PCA Analysis of Weston data
Data_file_PD= read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PD\\Weston PD Basline Reversal QC Mar 26 2018.csv',header=TRUE)
Data_file_PAL= read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\PAL\\Weston PAL Main Task Aggregated Mar 26 2018.csv',header=TRUE)
Data_file_5C= read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Mar 26 2018.csv',header=TRUE)

#PD

Data_PD = Data_file_PD[which((Data_file_PD$Task == 'PD Reversal 1') |(Data_file_PD$Task == 'PD Reversal 2') | (Data_file_PD$Task == 'PD Reversal 3')), ]
Data_PD = Data_PD[ ,c(3,4,5,6,7,8,9,11,18,17,49,81)]
Data_PD$No.Correction.Trials[is.na(Data_PD$No.Correction.Trials)] = 0
colnames(Data_PD) = c('Animal.ID','Test Site', 'Mouse.Strain','Genotype', 'Sex','Age','Task','Day','Accuracy','Corrections','Correct Latency','Reward Latency')
Data_PD =na.omit(Data_PD)
Data_PD$Age = as.factor(Data_PD$Age)
data.cols.pd <- Data_PD[c(10, 9, 11, 12)]
autoplot(prcomp(data.cols.pd), data = Data_PD, colour = 'Genotype', shape ='Sex', size ='Age')+ theme_classic()
res.pca.pd <- prcomp(data.cols.pd, scale = TRUE)
fviz_pca_var(res.pca.pd,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = ''
)


#PAL

Data_PAL = Data_file_PAL[, c(2,4,5,6,7,8,10,16,15,239,315)]
colnames(Data_PAL) = c('Animal.ID', 'Mouse.Strain','Genotype', 'Sex','Age','Task','Week','Accuracy','Corrections','Correct Latency','Reward Latency')
Data_PAL$Age = as.factor(Data_PAL$Age)
Data_PAL = na.omit(Data_PAL)
data.cols.pal <- Data_PAL[c(8,9,10,11)]
autoplot(prcomp(data.cols.pal), data = Data_PAL, colour = 'Genotype', shape ='Sex', size ='Age')+ theme_classic()

res.pca.pal <- prcomp(data.cols.pal, scale = TRUE)
fviz_pca_var(res.pca.pal,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = ''
)


#5C

Data_5C = Data_file_5C[ ,c(3,4,5,6,7,8,9,16,17,18,19,122,70)]
colnames(Data_5C) = c('Animal.ID','TestSite' ,'Mouse.Strain','Genotype', 'Sex','Age','StimulusDuration','Accuracy','Omissions','Premature','Perseverative','Correct Latency','Reward Latency')

Data_5C = na.omit(Data_5C)

data.cols.5c <- Data_5C[c( 8, 9, 10, 11, 12, 13)]
autoplot(prcomp(data.cols.5c), data = Data_5C, colour = 'Genotype', shape ='Sex', size ='Age')+ theme_classic()

res.pca.5c <- prcomp(data.cols.5c, scale = TRUE)
fviz_pca_var(res.pca.5c,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             title = ''
)
rm(list=ls(all=TRUE))
