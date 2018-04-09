## Library ##
library(reshape2)
library(car)
## Functions ##
Vigilance.ANOVA.Function = function(dataset){
  dataset = dataset[ ,c(3,4,5,6,7,8,9,14:23)]
  melt.data = melt(dataset, id.vars = c(1,2,3,4,5,6,7))
  new.cols = str_split_fixed(melt.data[ ,8], '.Block.',n=2)
  melt.data$variable = new.cols[ ,1]
  melt.data$Block = new.cols[ ,2]
  cast.data = dcast(melt.data, AnimalID + TestSite + Mouse.Strain + Genotype + Sex + Age.Months ~ variable + Stimulus.Length + Block,fun.aggregate = mean, na.rm=TRUE, value.var="value")
  
  idata = unique(melt.data[c('Stimulus.Length','Block')])
  idata = idata[order(idata$Stimulus.Length), ]
  idata$Stimulus.Length = as.factor(idata$Stimulus.Length)
  idata$Block = as.factor(idata$Block)
  
  strain.list = as.vector(unique(cast.data$Genotype))
  sex.list = as.vector(unique(cast.data$Sex))
  age.list = as.vector(unique(cast.data$Age.Months))
  
  final.datalist = list()
  for(a in 1:length(age.list)){
    for(b in 1:length(strain.list)){
      for(c in 1:length(sex.list)){
        analysis.data = cast.data[which((cast.data$Age.Months == age.list[a]) & (cast.data$Genotype == strain.list[b]) & (cast.data$Sex == sex.list[c])), ]
        accuracy.data = analysis.data[ ,c(1:6,7:26)]
        accuracy.matrix = as.matrix(accuracy.data[ ,c(7:26)])
        accuracy.lm = lm(accuracy.matrix ~ 1, data=accuracy.data)
        omission.data = analysis.data[ ,c(1:6,27:46)]
        omission.matrix = as.matrix(omission.data[ ,c(7:26)])
        omission.lm = lm(omission.matrix ~ 1, data=omission.data)
        accuracy.anova = Anova(accuracy.lm,idata=idata,idesign=~Stimulus.Length*Block,type='III')
        omission.anova = Anova(omission.lm,idata=idata,idesign=~Stimulus.Length*Block,type='III')
        final.datalist[[strain.list[b]]][[as.character(age.list[a])]][[sex.list[c]]][['Accuracy']]=accuracy.anova
        final.datalist[[strain.list[b]]][[as.character(age.list[a])]][[sex.list[c]]][['Omission']]=omission.anova
      }
    }
  }
  return(final.datalist)
}

## Read File ##
raw.file = file.choose()
raw.data = read.csv(raw.file)

## Prepare ANOVA ##
anova.data = Vigilance.ANOVA.Function(raw.data)

## Final Summary Table ##
strain.list = as.vector(names(anova.data))
measure.list = as.vector(names(anova.data$`5xFAD`$`4`$Female))
sex.list = as.vector(names(anova.data$`5xFAD`$`4`))
age.list = as.vector(names(anova.data$`5xFAD`))
template.file = summary(anova.data[[1]][[1]][[1]][[1]], multivariate=FALSE)
template.rownames = rownames(template.file[[4]])
template.rownames = template.rownames[2:length(template.rownames)]
template.pvalads = rownames(template.file[[5]])

hm.rownames = template.rownames
hm.rownames = gsub('Stimulus.Length','Stimulus Length',hm.rownames)
hm.rownames = gsub(':','*',hm.rownames)

strain.count = length(strain.list)
age.count = length(age.list)
measure.count = length(measure.list)
analysis.count = length(hm.rownames)
sex.count = length(sex.list)


summary.table = as.data.frame(matrix(nrow=(measure.count*analysis.count*strain.count*age.count*sex.count),ncol=10))
colnames(summary.table) = c('Strainr','Ager','Sexr','Measurer','Analysisr','df1r','df2r','Fr','pr','partial eta^2r')
row.modifier = 1
for(a in 1:strain.count){
  for(b in 1:age.count){
    for(c in 1:sex.count){
      for(d in 1:measure.count){
        temp.summary = summary(anova.data[[strain.list[a]]][[age.list[b]]][[sex.list[c]]][[measure.list[d]]], multivariate=FALSE)
        temp.main = temp.summary[[4]]
        temp.pvalad = temp.summary[[5]]
        for(e in 1:length(template.pvalads)){
          temp.main[which(rownames(temp.main) == template.pvalads[e]),6] = temp.pvalad[which(rownames(temp.pvalad) == template.pvalads[e]),2]
        }
        for(e in 1:length(hm.rownames)){
          partial.eta = temp.main[(e+1),1] / (temp.main[(e+1),1] + temp.main[(e+1),3])
          summary.table[row.modifier,1] = strain.list[a]
          summary.table[row.modifier,2] = age.list[b]
          summary.table[row.modifier,3] = sex.list[c]
          summary.table[row.modifier,4] = measure.list[d]
          summary.table[row.modifier,5] = hm.rownames[e]
          summary.table[row.modifier,6] = temp.main[(e+1),2]
          summary.table[row.modifier,7] = temp.main[(e+1),4]
          summary.table[row.modifier,8] = round(temp.main[(e+1),5],digits=2)
          summary.table[row.modifier,9] = round(temp.main[(e+1),6],digits=3)
          summary.table[row.modifier,10] = round(partial.eta, digits = 2)
          row.modifier = row.modifier + 1
        }
      }
    }
  }
}
