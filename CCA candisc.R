## Library ##
library(car)
library(dplyr)
library(candisc)
library(heplots)

## Functions ##

heplot.cancor <- function (
  mod,		         # output object from cancor
  which=1:2,       # canonical dimensions to plot
  scale,           # scale factor(s) for variable vectors in can space
  asp=1,           # aspect ratio, to ensure equal units?
  var.vectors = "Y", # which variable vectors to show?
  var.col=c("blue", "darkgreen"),  # colors for Y and X variable vectors and labels
  var.lwd=par("lwd"),
  var.cex=par("cex"),
  var.xpd=TRUE,     # was: par("xpd"),
  prefix = "Ycan",  # prefix for labels of canonical dimensions
  suffix = TRUE,   # add label suffix with can % ?
  terms=TRUE,  # terms to be plotted in canonical space / TRUE=all
  ...              # extra args passed to heplot
) {
  
  if (!inherits(mod, "cancor")) stop("Not a cancor object")
  if (mod$ndim < 2 || length(which)==1) {
    # using stop() here would terminate heplot.candiscList
    message("Can't do a 1 dimensional canonical HE plot")
    #	   plot(mod, which=which, var.col=var.col, var.lwd=var.lwd, prefix=prefix, suffix=suffix, ...) 
    return()
  }
  
  Yvars <- mod$names$Y
  scores <- data.frame(mod$scores$X, mod$scores$Y)
  scores <- data.frame(scores, mod$X)   # append X variables
  Xcoef <- mod$coef$X
  Ycoef <- mod$coef$Y
  Ycan <- colnames(Ycoef)
  
  canr <- mod$cancor
  lambda <- canr^2 / (1-canr^2)
  pct = 100*lambda / sum(lambda)
  
  if ((is.logical(terms) && terms) || terms=="X") {
    terms <- mod$names$X
  }
  # allow plotting the Xcan variables
  else if (length(terms)==1 && terms=="Xcan") terms=colnames(Xcoef)
  
  # make sure that all terms are available	
  if (!all(terms %in% colnames(scores))) {
    stop(paste(setdiff(terms, colnames(scores) ), "are not among the available variables"))
  }
  
  ##   Construct the model formula to fit mod$Yscores ~ Xscores in original lm()
  ##   using the mod$scores data.frame
  #browser()
  txt <- paste( "lm( cbind(",
                paste(Ycan, collapse = ","),
                ") ~ ",
                paste( terms, collapse = "+"), ", data=scores)" )
  can.mod <- eval(parse(text=txt))
  
  ##   Construct labels for canonical variables
  canvar <- Ycan[which]   # names of canonical variables to plot
  if (is.logical(suffix) & suffix)
    suffix <- paste( " (", round(pct[which],1), "%)", sep="" ) else suffix <- NULL
  canlab <- paste(prefix, which, suffix, sep="")
  
  ellipses <- heplot(can.mod, terms=terms, 
                     xlab=canlab[1], ylab=canlab[2], asp=asp, ...)
  
  struc <- mod$structure
  Xstructure <- struc$X.yscores[,which]
  Ystructure <- struc$Y.yscores[,which]
  vec <- rbind(
    if ("Y" %in% var.vectors) Ystructure else NULL,
    if ("X" %in% var.vectors) Xstructure else NULL)
  
  maxrms <- function(x) { max(sqrt(apply(x^2, 1, sum))) }
  if (missing(scale)) {
    vecrange <- range(vec)
    ellrange <- lapply(ellipses, range)
    vecmax <- maxrms(vec)
    ellmax <- max( maxrms(ellipses$E), unlist(lapply(ellipses$H, maxrms)) )
    scale <- floor(  0.9 * ellmax / vecmax )
    cat("Vector scale factor set to ", scale, "\n")
  }
  
  if ("Y" %in% var.vectors) vectors(Ystructure, scale=scale, col=var.col[1], lwd=var.lwd, cex=var.cex, xpd=var.xpd)
  if ("X" %in% var.vectors) vectors(Xstructure, scale=scale, col=var.col[2], lwd=var.lwd, cex=var.cex, xpd=var.xpd)
  
  invisible(ellipses)
  
  
}

## Input Raw Data ##
raw.data.probe = read.csv('C:\\Users\\dpalmer\\Documents\\Weston_R_Script\\Data\\Raw\\5CSRTT\\Weston 5CSRTT Probe Aggregated QC Oct 12 2017 Updated.csv')

## Rescore Text vars ##

raw.data.probe$TestSite = dplyr::recode(raw.data.probe$TestSite,'Site1'=1,'Site2'=2)

raw.data.probe$Genotype = dplyr::recode(raw.data.probe$Genotype,"C57BL6"=1,"B6SJLF1/J"=1,"B6129SF2/J"=1,"APPPS1"=2,"5xFAD"=2,"3xTG-AD"=2)

raw.data.probe$Sex = dplyr::recode(raw.data.probe$Sex, 'Female'=1,'Male'=2)

## Generate Two Groupings for Canonical Correlation ##
strain.data = list()
strain.data$APP = raw.data.probe[which(raw.data.probe[ ,5] == 'APP-PS1'),c(3:9,16:19,70,122)]
colnames(strain.data$APP) = c('Animal.ID','TestSite','Mouse Strain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature Responses','Perseverative Responses','Reward Collection Latency','Correct Response Latency')
strain.data$TG5x = raw.data.probe[which(raw.data.probe[ ,5] == '5xFAD'),c(3:9,16:19,70,122)]
colnames(strain.data$TG5x) = c('Animal.ID','TestSite','Mouse Strain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature Responses','Perseverative Responses','Reward Collection Latency','Correct Response Latency')
strain.data$TG3x = raw.data.probe[which(raw.data.probe[ ,5] == '3xTG-AD'),c(3:9,16:19,70,122)]
colnames(strain.data$TG3x) = c('Animal.ID','TestSite','Mouse Strain','Genotype','Sex','Age','StimulusLength','Accuracy','Omissions','Premature Responses','Perseverative Responses','Reward Collection Latency','Correct Response Latency')

## Generate List for Canonical Factors
canon.raw.list = list()
canon.raw.list$APP$MouseFactor = strain.data$APP[ ,c(2,4,5,6,7)]
canon.raw.list$APP$BehaviourMeasure = strain.data$APP[ ,c(8,9,10,11,12,13)]
canon.raw.list$TG5x$MouseFactor = strain.data$TG5x[ ,c(2,4,5,6,7)]
canon.raw.list$TG5x$BehaviourMeasure = strain.data$TG5x[ ,c(8,9,10,11,12,13)]
canon.raw.list$TG3x$MouseFactor = strain.data$TG3x[ ,c(2,4,5,6,7)]
canon.raw.list$TG3x$BehaviourMeasure = strain.data$TG3x[ ,c(8,9,10,11,12,13)]

## CCA - candisc ##

cca.list = list()
cca.list$APP = cancor(canon.raw.list$APP$MouseFactor,canon.raw.list$APP$BehaviourMeasure,set.names=c('Mouse Factors','Behavioural Measures'))
cca.list$TG5x = cancor(canon.raw.list$TG5x$MouseFactor,canon.raw.list$TG5x$BehaviourMeasure,set.names=c('Mouse Factors','Behavioural Measures'))
cca.list$TG3x = cancor(canon.raw.list$TG3x$MouseFactor,canon.raw.list$TG3x$BehaviourMeasure,set.names=c('Mouse Factors','Behavioural Measures'))
