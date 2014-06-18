
source("load.R")
library(dplyr)
library(sciplot)


# BH pipemodel relationships
datl <- studyWithVars(baad, c("a.stbh", "a.lf","m.lf"), returnwhat="list")
dat <- rbind_all(datl)

# Roth2007 massive outlier, remove here
dat <- subset(dat, dataset != "Roth2007")

# remove DG; there is only one species which makes it confusing
dat <- subset(dat, pft != "DG")

pipe_bh <- summarize(group_by(dat, dataset, species), 
           alf_stbh=mean(a.lf / a.stbh),
           mlf_stbh=mean(m.lf / a.stbh),
           pft=first(pft),
           MAP=mean(MAP),
           MAT=mean(MAT),
           growingCondition=first(growingCondition),
           vegetation=first(vegetation))
           
# Basal pipemodel relationships
datl <- studyWithVars(baad, c("a.stba", "a.lf","m.lf"), returnwhat="list")
dat <- rbind_all(datl)

pipe_ba <- summarize(group_by(dat, dataset, species), 
                     alf_stba=mean(a.lf / a.stba),
                     mlf_stba=mean(m.lf / a.stba),
                     pft=first(pft),
                     MAP=mean(MAP),
                     MAT=mean(MAT),
                     growingCondition=first(growingCondition),
                     vegetation=first(vegetation))



BARplot <- function(x,y,ymult=1.2,...){
  
  N <- table(x)
  tap <- tapply(y,x,mean,na.rm=T)
  Ylim <- c(0, max(tap)*ymult)
  b <- bargraph.CI(x,y,ylim=Ylim,err.width=0.04,...)
  
  d <- 0.3
  text(c(0.3, b$xvals+d), max(tap)*(1.1/1.2)*ymult, labels=c("N =", as.character(unname(N))))
  
  return(invisible(b))
}

Cols <- c("blue","red","darkorange")

# Cols <- c("blue","royalblue","red","darkorange")
windows(8,8)
par(mfrow=c(2,2), mar=c(5,5,2,2))
with(pipe_bh, BARplot(pft,alf_stbh, xlab="Plant functional type", 
                       ylab=expression(Leaf~area/breast~height~stem~area~~(m^2~m^-2)),
                       col=Cols))
with(pipe_bh, BARplot(pft,mlf_stbh, xlab="Plant functional type",
                       ylab=expression(Leaf~mass/breast~height~stem~area~~(kg~m^-2)),
                       col=Cols
                       ))
with(pipe_ba, BARplot(pft,alf_stba, xlab="Plant functional type", 
                          ylab=expression(Leaf~area/basal~stem~area~~(m^2~m^-2)),
                          col=Cols))
b <- with(pipe_ba, BARplot(pft,mlf_stba, xlab="Plant functional type",
                          ylab=expression(Leaf~mass/basal~stem~area~~(kg~m^-2)),
                          col=Cols,ymult=1.3
))





#------------------------------------------------------------------------------------------#
plotit <- function(i){
  
  dataset <- datl[[i]]
  
  if(!("grouping" %in% names(datl[[i]])) || length(unique(datl[[i]]$grouping))==1){
    dataset$grouping <- 1
    palette(c("black","black"))
  } else {
    dataset$grouping <- as.factor(dataset$grouping)
    palette(rainbow(nlevels(dataset$grouping)))
  }
  
  par(mfrow=c(1,2))
  with(dataset, plot(log10(a.stbh), log10(m.lf), pch=19, col=grouping, axes=F))
  magaxis(side=1:2, unlog=1:2, box=T)
  with(dataset, plot(log10(a.stbh), log10(a.lf), pch=19, col=grouping, axes=F))
  magaxis(side=1:2, unlog=1:2, box=T)
  
  title(names(datl)[i], outer=T, line=-2, cex=0.9)
}


pdf("output/pdf/pipeplotbystudygroup.pdf", width=8, height=4)
for(i in 1:length(datl))plotit(i)
dev.off()




