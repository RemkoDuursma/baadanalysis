
source("R/prepareDataset.R")

plotPredBarplot <- function(model, ylab="", Hs=c(2,15), nboot=100,...){
  
  
  newdat <- with(dataset, expand.grid(pft=unique(dataset$pft), 
                                      lh.t=log10(Hs)))
  P <- bootMer(model, nsim=nboot, FUN=function(.)predict(., newdat, re.form=NA))
  sims <-  10^as.data.frame(P)
  
  newdat$ypred <- 10^predict(model, newdat, re.form=NA)
  
  m <- matrix(newdat$ypred, ncol=3, byrow=TRUE)
  colnames(m) <- unique(dataset$pft)
  
  qu <- apply(sims, 2, quantile, probs=c(0.025,0.975))
  
  b <- barplot(m, beside=TRUE, ylim=c(0,max(qu[2,])),ylab=ylab,...)
  f <- function(x)as.vector(matrix(x,ncol=3,byrow=T))
  arrows(x0=as.vector(b),
         y0=f(qu[1,]),
         y1=f(qu[2,]),
         angle=90,code=3,length=0.1)
}


fit5 <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
plotPredBarplot(fit5, ylab="m.lf / a.stbh")

fit6 <- lmer(lalf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
plotPredBarplot(fit6, ylab="a.lf / a.stbh")

fit8 <- lmer(lmlf_mso ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
plotPredBarplot(fit8, ylab="m.lf / m.so")

fit9 <- lmer(lalf_mso ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
plotPredBarplot(fit9, ylab="a.lf / m.so")

fit10 <- lmer(lmrt_mso ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
plotPredBarplot(fit10, ylab="m.rt / m.so")



#-------------------------------------------------------------------------------------------#

plotPredBarplot2 <- function(model, ylab="", Hs=c(2,15), nboot=100, dfr=dataset2, ...){
  
  
  newdat <- with(dfr, expand.grid(pft=unique(dfr$pft), 
                                      bortemptrop=unique(dfr$bortemptrop),
                                      lh.t=log10(Hs)
                                      ))
  
  P <- bootMer(model, nsim=nboot, FUN=function(.)predict(., newdat, re.form=NA))
  sims <-  10^as.data.frame(P)
  
  newdat$ypred <- 10^predict(model, newdat, re.form=NA)
  
  
  givenames <- function(x){
    colnames(x) <- levels(dfr$pft)
    rownames(x) <- rep(levels(dfr$bortemptrop),2)
    return(x)
  }
  
  m <- matrix(newdat$ypred, ncol=3, byrow=T)
  m <- givenames(m)
  
  qu <- apply(sims, 2, quantile, probs=c(0.025,0.975))
  qul <- matrix(qu[1,],ncol=3,byrow=T)
  quu <- matrix(qu[2,],ncol=3,byrow=T)
  qul <- givenames(qul)
  quu <- givenames(quu)
  
  # Remove combo's not in dataset.
  setCellsToNA <- function(x, cells=rbind(c(1,3,4,6),c(2,3,2,3))){
    for(i in 1:ncol(ii))x[ii[1,i],ii[2,i]] <- NA
    return(x)
  }
  quu <- setCellsToNA(quu)
  qul <- setCellsToNA(qul)
  m <- setCellsToNA(m)
  
  b <- barplot2(m, beside=TRUE, plot.ci=TRUE, ci.l=qul, ci.u=quu, ci.width=0.15,
           col=c(rep("lightgrey",3), rep("dimgrey",3)))
  
  b <- as.vector(b)
  x <- rep(c("Bo","Te","Tr"), 6)
  mtext(side=1, line=-0.2, text=x, cex=0.8, at=b)
  legend("topleft", c(paste("H =",Hs)), fill=c("lightgrey","dimgrey"), cex=0.6, bty='n')

  
  newdat$ypred[newdat$pft == "EG" & newdat$bortemptrop == "tropical"] <- NA
  newdat$ypred[newdat$pft == "EA" & newdat$bortemptrop == "boreal"] <- NA
  
return(invisible(newdat))
}


fita <- lmer(lmlf_astbh ~ pft + bortemptrop + lh.t + pft:lh.t + bortemptrop:lh.t + (lh.t|Group), 
             data=dataset2)
fitb <- lmer(lalf_astbh ~ pft + bortemptrop + lh.t + pft:lh.t + bortemptrop:lh.t + (lh.t|Group), 
             data=dataset2)

fitc <- lmer(lsla ~ pft + bortemptrop + lh.t + pft:lh.t + bortemptrop:lh.t + (lh.t|Group), 
             data=dataset2)

# fita.0 <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), 
#              data=dataset2)
# anova(fita.0, fita)


data <- plotPredBarplot2(fita, Hs=c(2,10), nboot=10)
datb <- plotPredBarplot2(fitb, Hs=c(2,10), nboot=10)
datc <- plotPredBarplot2(fitc, Hs=c(2,10), nboot=10)

names(data)[4] <- "mlf_astbh"
names(datb)[4] <- "alf_astbh"
names(datc)[4] <- "sla" 
dat <- merge(data, datb)
dat <- merge(dat, datc)

Cols <- c("blue","red","darkorange")
palette(Cols)

datsmall <- subset(dat, lh.t  == min(lh.t))
datlarge <- subset(dat, lh.t  == max(lh.t))

par(xaxs="i", yaxs="i")
pchs <- c(15,17,19)
with(datsmall, plot(1/sla, mlf_astbh, ylim=c(0,1000),xlim=c(0,0.4),
                    xlab=expression(LMA~(kg~m^-2)),
                    ylab=expression(W[F]/A[S]~~(kg~m^-2)),
               col=pft, pch=pchs[bortemptrop]))
abline(lm(mlf_astbh ~ I(1/sla), data=datsmall))
pchs <- c(22,24,21)
with(datlarge, points(1/sla, mlf_astbh, 
                    col=pft, pch=pchs[bortemptrop]))
abline(lm(mlf_astbh ~ I(1/sla), data=datlarge))
legend("topleft", levels(dat$bortemptrop), pch=pchs, cex=0.7, pt.cex=1, bty='n')


par(xaxs="i", yaxs="i")
pchs <- c(15,17,19)
with(datsmall, plot(1/sla, alf_astbh, xlim=c(0,0.4),ylim=c(0,5000),
                    xlab=expression(LMA~(kg~m^-2)),
                    ylab=expression(A[L]/A[S]~~(kg~m^-2)),
                    col=pft, pch=pchs[bortemptrop]))
# abline(lm(alf_astbh ~ I(1/sla), data=datsmall))
pchs <- c(22,24,21)
with(datlarge, points(1/sla, alf_astbh, 
                      col=pft, pch=pchs[bortemptrop]))
# abline(lm(alf_astbh ~ I(1/sla), data=datlarge))
legend("topright", levels(dat$bortemptrop), pch=pchs, cex=0.7, pt.cex=1, bty='n')






