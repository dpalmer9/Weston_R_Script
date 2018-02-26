## Prepare ANOVA Table ##
strain.list = as.vector(names(probe.anova))
measure.list = as.vector(names(probe.anova$APP$Month4))
measure.list = measure.list[c(3:6)]
age.list = as.vector(names(probe.anova$APP))
template.file = summary(probe.anova[[1]][[1]][[1]], multivariate=FALSE)
template.rownames = rownames(template.file[[4]])
template.rownames = template.rownames[2:length(template.rownames)]
template.pvalads = rownames(template.file[[5]])

hm.rownames = template.rownames
hm.rownames = gsub('Genotype','G',hm.rownames)
hm.rownames = gsub('Sex','Sx',hm.rownames)
hm.rownames = gsub('Week','B',hm.rownames)
hm.rownames = gsub(':','*',hm.rownames)

strain.count = length(strain.list)
age.count = length(age.list)
measure.count = length(measure.list)
analysis.count = length(hm.rownames)


map.list = list()
summary.table = as.data.frame(matrix(nrow=(strain.count*age.count*measure.count*analysis.count),ncol=9))
colnames(summary.table) = c('Strain','Age','Measure','Analysis','df1','df2','F','p','partial eta^2')
for(a in 1:length(strain.list)){
  for(b in 1:length(age.list)){
    for(d in 1:length(measure.list)){
      temp.summary = summary(probe.anova[[strain.list[a]]][[age.list[b]]][[measure.list[d]]], multivariate=FALSE)
      temp.main = temp.summary[[4]]
      temp.pvalad = temp.summary[[5]]
      for(c in 1:length(template.pvalads)){
        temp.main[which(rownames(temp.main) == template.pvalads[c]),6] = temp.pvalad[which(rownames(temp.pvalad) == template.pvalads[c]),2]
      }
      for(c in 1:length(hm.rownames)){
        row.modifier = (((a-1)*age.count*measure.count*analysis.count) + ((b-1)*measure.count*analysis.count) + ((d-1)*analysis.count) + c)
        partial.eta = temp.main[(c+1),1] / (temp.main[(c+1),1] + temp.main[(c+1),3])
        summary.table[row.modifier,1] = strain.list[a]
        summary.table[row.modifier,2] = age.list[b]
        summary.table[row.modifier,3] = measure.list[d]
        summary.table[row.modifier,4] = hm.rownames[c]
        summary.table[row.modifier,5] = temp.main[(c+1),2]
        summary.table[row.modifier,6] = temp.main[(c+1),4]
        summary.table[row.modifier,7] = round(temp.main[(c+1),5],digits=2)
        summary.table[row.modifier,8] = round(temp.main[(c+1),6],digits=3)
        summary.table[row.modifier,9] = round(partial.eta, digits = 2)
      }
    }
  }
}
map.list[['Summary']] = summary.table