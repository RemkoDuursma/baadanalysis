
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

plotPredBarplot2 <- function(model, ylab="", Hs=c(2,15), nboot=100,...){
  
  
  newdat <- with(dataset, expand.grid(pft=unique(dataset2$pft), 
                                      bortemptrop=unique(dataset2$bortemptrop),
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


fita <- lmer(lmlf_astbh ~ pft + bortemptrop + lh.t + pft:lh.t + bortemptrop:lh.t + (1|Group), data=dataset2)


