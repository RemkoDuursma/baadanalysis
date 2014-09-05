
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



#--------------------------------------------------------------------------------------#

# pftlong = pft * bortemptrop
plotPredBarplot3 <- function(model, yvar="ypred", ylab="", Hs=c(2,15), nboot=100, dfr=dataset2,...){
  
  n <- length(unique(dfr$pftlong))
  newdat <- with(dfr, expand.grid(pftlong=unique(dfr$pftlong), 
                                      lh.t=log10(Hs)))
  P <- bootMer(model, nsim=nboot, FUN=function(.)predict(., newdat, re.form=NA))
  sims <-  10^as.data.frame(P)
  
  newdat$ypred <- 10^predict(model, newdat, re.form=NA)
  
  m <- matrix(newdat$ypred, ncol=n, byrow=TRUE)
  colnames(m) <- unique(dfr$pftlong)
  
  qu <- apply(sims, 2, quantile, probs=c(0.025,0.975))
  
  b <- barplot(m, beside=TRUE, ylim=c(0,max(qu[2,])),ylab=ylab,...)
  f <- function(x)as.vector(matrix(x,ncol=n,byrow=T))
  arrows(x0=as.vector(b),
         y0=f(qu[1,]),
         y1=f(qu[2,]),
         angle=90,code=3,length=0.1)
  
  newdat$lci <- qu[1,]
  newdat$uci <- qu[2,]
  names(newdat)[3:5] <- c(yvar, paste0(yvar,"L"), paste0(yvar,"U"))
  return(invisible(newdat))
}


fita <- lmer(lmlf_astbh ~ pftlong + lh.t + pftlong:lh.t + (lh.t|Group), 
             data=dataset2)

fitb <- lmer(lalf_astbh ~ pftlong + lh.t + pftlong:lh.t + (lh.t|Group), 
             data=dataset2)

fitc <- lmer(lsla ~ pftlong + lh.t + pftlong:lh.t + (lh.t|Group), 
             data=dataset2)

fitd <- lmer(lmlf_mso ~ pftlong + lh.t + pftlong:lh.t + (lh.t|Group), 
             data=dataset2)

fite <- lmer(lalf_mso ~ pftlong + lh.t + pftlong:lh.t + (lh.t|Group), 
             data=dataset2)

N <- 10
data <- plotPredBarplot3(fita, "lmlf_astbh", Hs=c(2,10), nboot=N)

datb <- plotPredBarplot3(fitb, "lalf_astbh", Hs=c(2,10), nboot=N, dfr=droplevels(subset(dataset2, !is.na(lalf_astbh) & !is.na(lh.t))))
datc <- plotPredBarplot3(fitc, "lsla", Hs=c(2,10), nboot=N,
                         dfr=droplevels(subset(dataset2, !is.na(lsla) & !is.na(lh.t))))
datd <- plotPredBarplot3(fitd, "lmlf_mso", Hs=c(2,10), nboot=N,
                         dfr=droplevels(subset(dataset2, !is.na(lmlf_mso) & !is.na(lh.t))))
date <- plotPredBarplot3(fite, "lalf_mso", Hs=c(2,10), nboot=N,
                         dfr=droplevels(subset(dataset2, !is.na(lalf_mso) & !is.na(lh.t))))

dat <- merge(data, datb, all=T)
dat <- merge(dat, datc, all=T)
dat <- merge(dat, datd, all=T)
dat <- merge(dat, date, all=T)

dat$lma <- 1/dat$lsla
dat$lmaU <- 1/dat$lslaU
dat$lmaL <- 1/dat$lslaL

dat1 <- subset(dat, lh.t == min(lh.t))
dat2 <- subset(dat, lh.t == max(lh.t))



pc <- function(xvar, yvar,dfr){

    arrows(x0=dfr[,paste0(xvar,"L")],
           x1=dfr[,paste0(xvar,"U")],
         y0=dfr[,yvar],
         y1=dfr[,yvar],
         angle=90,code=3,length=0.05)
    arrows(x0=dfr[,xvar],
           x1=dfr[,xvar],
           y0=dfr[,paste0(yvar,"L")],
           y1=dfr[,paste0(yvar,"U")],
           angle=90,code=3,length=0.05)
  
}


Cols <- c("dimgrey", "royalblue")

par(mfrow=c(1,2),xaxs="i", yaxs="i")
with(dat1, plot(lma, lmlf_astbh, ylim=c(0,1000), xlim=c(0,0.4), pch=19, col=Cols[1],
panel.first={abline(lm(lmlf_astbh ~ lma-1, data=dat1),col="grey")
             abline(lm(lmlf_astbh ~ lma-1, data=dat2),col="grey")
            pc("lma","lmlf_astbh",dat1)
             pc("lma","lmlf_astbh",dat2)}))
with(dat2, points(lma, lmlf_astbh, pch=19, col=Cols[2]))

with(dat1, plot(1/lsla, lalf_astbh, xlim=c(0,0.4), ylim=c(0,5000), pch=19, col=Cols[1],
  panel.first={pc("lsla","lalf_astbh",dat1)
               pc("lsla","lalf_astbh",dat2)}))
with(dat2, points(1/lsla, lalf_astbh, pch=19, col=Cols[2]))



par(mfrow=c(1,2),xaxs="i", yaxs="i")
with(dat1, plot(lma, lmlf_mso, xlim=c(0,0.4), ylim=c(0,0.4), pch=19, col=Cols[1],
  panel.first={abline(lm(lmlf_mso ~ lma-1, data=dat1),col="grey")
               abline(lm(lmlf_mso ~ lma-1, data=dat2),col="grey")
               pc("lma","lmlf_mso",dat1)
               pc("lma","lmlf_mso",dat2)}))
with(dat2, points(lma, lmlf_mso, pch=19, col=Cols[2]))


with(dat1, plot(lma, lalf_mso, ylim=c(0,4), xlim=c(0,0.4), pch=19, col=Cols[1],
                panel.first={pc("lma","lalf_mso",dat1)
                             pc("lma","lalf_mso",dat2)}))
with(dat2, points(lma, lalf_mso, pch=19, col=Cols[2]))







with(dataset2, plot(lh.t, lmlf_mso))










