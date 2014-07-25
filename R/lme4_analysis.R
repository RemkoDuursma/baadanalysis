

source("load.R")

library(lme4)
library(lmerTest)
library(visreg)
library(LMERConvenienceFunctions)

dataset <- baad
dataset$Group <- paste(dataset$studyName, dataset$speciesMatched)
dataset <- droplevels(subset(dataset, growingCondition %in% c("FW","PM","PU","FE") & pft != "DG"))
dataset$pft <- as.factor(dataset$pft)

pftcols <- list(EA="red", EG="hotpink", DA="blue")




#-----------------------------------------------------------------------------#

dataset <- subset(dataset, !is.na(MAT) & !is.na(pft))

dataset <- subset(dataset, studyName != "Roth2007")

# Nothing
fit0 <- lmer(log10(m.lf) ~ log10(a.stbh) + (log10(a.stbh)|Group),
             data=dataset)

# Add pft
fit1 <- lmer(log10(m.lf) ~ pft*log10(a.stbh) + (log10(a.stbh)|Group),
data=dataset)

library(visreg)
visreg(fit1, "a.stbh", by="pft", xtrans=log10, overlay=T, partial=FALSE)


# Add vegetation
fit2 <- lmer(log10(m.lf) ~ vegetation + pft + pft:log10(a.stbh) + (log10(a.stbh)|Group),
             data=dataset)

visreg(fit2, "a.stbh", by="vegetation", xtrans=log10, overlay=T, partial=FALSE)


# 
fit3 <- lmer(log10(m.lf) ~ MAT + pft + pft:log10(a.stbh) + (log10(a.stbh)|Group),
             data=dataset)

visreg(fit2, "a.stbh", by="vegetation", xtrans=log10, overlay=T, partial=FALSE)


fit4 <- lmer(log10(m.lf) ~ MAT*log10(a.stbh) + (log10(a.stbh)|Group),
             data=dataset)


dataset$lmlf_astbh <- with(dataset, log10(m.lf/a.stbh))
dataset$lalf_astbh <- with(dataset, log10(a.lf/a.stbh))
dataset$lh.t <- with(dataset, log10(h.t))

fit5 <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)
fit6 <- lmer(lalf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)

visreg(fit5, "lh.t", by="pft", overlay=T)
visreg(fit6, "lh.t", by="pft", overlay=T)


fit5 <- lmer(lmlf_astbh ~ pft + lh.t + pft:lh.t + (1|Group), data=dataset)

#-----------------------------------------------------------------------------#

# prediction barplots.
# because log-log plots show so little, which is why others have jumped to
# universal scaling etc. PFT has a large effect but it is so hard to see on a plot where
# ht varies from 0.01 to 100m



plotPredBarplot <- function(model, ylab="", Hs=c(2,15), nboot=100){

  newdat <- with(dataset, expand.grid(pft=unique(dataset$pft), 
                                      lh.t=log10(Hs)))
  P <- bootMer(model, nsim=nboot, FUN=function(.)predict(., newdat, re.form=NA))
  sims <-  10^as.data.frame(P)
  
  newdat$ypred <- 10^predict(model, newdat, re.form=NA)
  
  m <- matrix(newdat$ypred, ncol=3, byrow=TRUE)
  colnames(m) <- unique(dataset$pft)
  
  qu <- apply(sims, 2, quantile, probs=c(0.025,0.975))
  
  b <- barplot(m, beside=TRUE, ylim=c(0,max(qu[2,])))
  f <- function(x)as.vector(matrix(x,ncol=3,byrow=T))
  arrows(x0=as.vector(b),
         y0=f(qu[1,]),
         y1=f(qu[2,]),
         angle=90,code=3,length=0.1)
}


plotPredBarplot(fit5, ylab="m.lf / a.stbh")
plotPredBarplot(fit6, ylab="a.lf / a.stbh")


#-----------------------------------------------------------------------------#
runit <- function(xvar,yvar){

  fitit <- function(xvar, yvar){
  
    dataset$Y <<- log10(dataset[,yvar])
    dataset$X <<- log10(dataset[,xvar])
    fito <- lmer(Y ~ pft*X + (X|Group), data=dataset)
    
    return(fito)
  }
  
  p <- fitit(xvar,yvar)
  
  with(dataset, plot(X,Y, pch=16, col="grey", cex=0.8,
                     xlab=xvar, ylab=yvar, axes=FALSE))
  magaxis(side=1:2, unlog=1:2)
  m <- matrix(fixef(p), ncol=3, byrow=T)
  colnames(m) <- levels(dataset$pft)
  for(i in 2:ncol(m))m[,i] <- m[,i] + m[,1]
  for(i in 1:ncol(m)){
    b0 <- m[1,i]
    b1 <- m[2,i]
    abline(b0,b1,col=pftcols[[levels(dataset$pft)[i]]])
  }
  return(invisible(list(model=p, coef=m)))
}



m1 <- runit("h.t","m.st")
m2 <- runit("a.stbh","h.t")

dataset$acs <- with(dataset, d.cr * c.d)
m3 <- runit("h.t", "acs")

m4 <- runit("acs", "a.lf")

makedf <- function(x, varnames){
  c1 <- coef(x$model)$Group[,c("(Intercept)","X")]
  names(c1) <- varnames
  c1$Group <- rownames(c1)
return(c1)
}

dfr <- merge(makedf(m1,c("c1","c2")), makedf(m2, c("c3","c4")))
dfr <- merge(dfr, dataset[,c("Group","pft")], all=F)
dfr <- dfr[!duplicated(dfr),]





fitit <- function(xvar, yvar){
  
  dataset$Y <- log10(dataset[,yvar])
  dataset$X <- log10(dataset[,xvar])
  fito <- lmer(Y ~ pft*X + (X|Group), data=dataset)
  
  return(fito)
}
fitit2 <- function(xvar, yvar){
  
  dataset$Y <- log10(dataset[,yvar])
  dataset$X <- log10(dataset[,xvar])
  fito <- lme(Y ~ pft*X, random=~X|Group, data=dataset, na.action=na.omit)
  
  return(fito)
}


f1 <- fitit("a.stbh","h.t")
f2 <- fitit2("a.stbh","h.t")

