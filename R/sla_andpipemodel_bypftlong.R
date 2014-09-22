
source("load.R")
source("R/preparedataset.R")


meansbypftlong <- function(yvar1, yvar2, pftvar, xlab=expression(Specific~leaf~area~(m^2~kg^-1)), 
                           xlim=c(0,25),main="",
                           ylab1=NULL, ylab2=NULL, addtrend=c(FALSE,FALSE),
                           ylim1=NULL, ylim2=NULL, dataset){
  
  if(pftvar == "pftlong"){
    dat <- droplevels(subset(dataset, pftlong %in% 
                             c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                               "EG-temperate")))
  } else {
    dat <- dataset
  }
  
  dat$Y1 <- dat[,yvar1]
  dat$Y2 <- dat[,yvar2]
  dat$P <- dat[,pftvar]
  
  mixmean <- function(yvar, dataset=dat){
    
    dataset$Y <- dataset[,yvar]
    f <- lmer(Y ~ P - 1 + (1|Group), data=dataset)
    
    ci <- suppressMessages(confint(f)[-c(1,2),])
    ci <- as.data.frame(cbind(10^fixef(f),10^ci))
    rownames(ci) <- gsub("P","", rownames(ci))
    names(ci) <- c("y","lci","uci")
    return(ci)
  }
  
  sla <- mixmean("lsla")
  mlfastbh <- mixmean(yvar1)
  alfastbh <- mixmean(yvar2)
  
  palette(rainbow(length(unique(dat$P))))
  
  # SHOULD be updated
  lmfit1 <- lm(mlfastbh$y ~ sla$y)
  lmfit2 <- lm(alfastbh$y ~ sla$y)
  
  # model fit for multiple comparisons
  f1 <- lmer(Y1 ~ P - 1 + (1|Group), data=dat)
  f2 <- lmer(Y2 ~ P - 1 + (1|Group), data=dat)
  
  # Same results, much faster.
  lets1 <- cld(glht(f1, linfct=mcp(P="Tukey")))$mcletters$Letters
  lets2 <- cld(glht(f2, linfct=mcp(P="Tukey")))$mcletters$Letters

  #lets1 <- lets1[order(sla$y)]
  #lets2 <- lets2[order(sla$y)]
  
  #windows(4,7)
  par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
  plot(sla$y, mlfastbh$y, xlim=xlim,axes=FALSE, pch=19, col=1:5, cex=1.3,
       ylim=ylim1,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=mlfastbh$y, y1=mlfastbh$y,code=3,angle=90,length=0.025,col=1:5)
         arrows(x0=sla$y, x1=sla$y, y0=mlfastbh$lci, y1=mlfastbh$uci,code=3,angle=90,length=0.025,col=1:5)
         if(addtrend[1])ablinepiece(lmfit1)
       })
  axis(2)
  box()
  u <- par()$usr
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets1, pos=3, cex=0.9)
  axis(1,labels=FALSE)
  legend("topright", rownames(sla),pch=19,col=palette(), cex=0.7,pt.cex=1.2)
  plot(sla$y, alfastbh$y, xlim=xlim,pch=19, col=1:5, cex=1.3,
       ylim=ylim2,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=alfastbh$y, y1=alfastbh$y,code=3,angle=90,length=0.025,col=1:5)
         arrows(x0=sla$y, x1=sla$y, y0=alfastbh$lci, y1=alfastbh$uci,code=3,angle=90,length=0.025,col=1:5)
         if(addtrend[2])ablinepiece(lmfit2)
       })
  u <- par()$usr
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets2, pos=3, cex=0.9)
  mtext(side=1, text=xlab, line=3, outer=T)
  mtext(side=2, at = 0.25, text=ylab2, line=3, outer=T)
  mtext(side=2, at = 0.75, text=ylab1, line=3, outer=T)

  mtext(side=4, text=main, line=2, outer=TRUE)
  
}
  

# Pipe model, breast height
to.pdf(
  meansbypftlong("lmlf_astbh","lalf_astbh", "pftlong", ylab2=expression(A[L]/A[Sbh]~~(m^2~m^-2)),
    ylab1=expression(M[L]/A[Sbh]~~(kg~m^-2)), dataset=subset(dataset2, h.t > 1.3),
    ylim1=c(0,500),ylim2=c(0,4500)),width=4, height=7,
  filename="output/figures/piperatios_bypftlong_bh.pdf")

# Pipe model, basal height
to.pdf(
meansbypftlong("lmlf_astba","lalf_astba", ylab2=expression(A[L]/A[Sba]~~(m^2~m^-2)),
               ylim1=c(0,400), ylim2=c(0,3000),dataset=dataset2, #subset(dataset2, h.t > 1.3),
               ylab1=expression(M[L]/A[Sba]~~(kg~m^-2))),width=4, height=7,
  filename="output/figures/piperatios_bypftlong_ba.pdf")


# Pipe model, estimated and measured basal height
to.pdf(
  meansbypftlong("lmlf_astba2","lalf_astba2", ylab2=expression(A[L]/A[Sba~est]~~(m^2~m^-2)),
                 ylim1=c(0,400), ylim2=c(0,3000),dataset=subset(dataset2, h.t > 1.3),
                 ylab1=expression(M[L]/A[Sba~est]~~(kg~m^-2))),width=4, height=7,
  filename="output/figures/piperatios_bypftlong_ba_est.pdf")


# Basal diameter, small plants only
to.pdf(
  meansbypftlong("lmlf_astba","lalf_astba", ylab2=expression(A[L]/A[Sba]~~(m^2~m^-2)),
                 ylim1=c(0,400), ylim2=c(0,3000),xlim=c(0,40),
                 dataset=subset(dataset2, h.t < 1.3),
                 ylab1=expression(M[L]/A[Sba]~~(kg~m^-2))),width=4, height=7,
  filename="output/figures/piperatios_bypftlong_ba_smallonly.pdf")


# Root mass fraction, leaf mass fraction
to.pdf(
meansbypftlong("lmrt_mso","lmlf_mso", addtrend=c(F,F),dataset=subset(dataset2, h.t > 1.3),
               ylab1=expression(M[R]/(M[F] + M[W])~~(kg~kg^-1)),ylim1=c(0,0.5),
               ylab2=expression(M[F]/M[W]~~(kg~kg^-1)),ylim2=c(0,0.25)),width=4, height=7,
  filename="output/figures/RMF_LMF_bypftlong.pdf")





to.pdf({
meansbypftlong("lmlf_astbh","lalf_astbh", "pft", ylab2=expression(A[L]/A[Sbh]~~(m^2~m^-2)),
               ylab1=expression(M[L]/A[Sbh]~~(kg~m^-2)), dataset=subset(dataset2, h.t > 1.3),
               ylim1=c(0,500),ylim2=c(0,4000))

meansbypftlong("lmlf_astba","lalf_astba", "pft", ylab2=expression(A[L]/A[Sba]~~(m^2~m^-2)),
               ylab1=expression(M[L]/A[Sbh]~~(kg~m^-2)), dataset=subset(dataset2, h.t > 1.3),
               ylim1=c(0,500),ylim2=c(0,4000))
}





# Histograms with pipe model ratios.
dat <- subset(dataset2, h.t > 1.3)
d <- split(dat, dat$pftlong)



histbypftlong <- function(yvar, dataset, nbin=100, log=TRUE, col=1:5,
                          xlab=NULL, ylab="Nr. individuals", meanline=TRUE, openwin=FALSE){
  
  
  yall <- eval(substitute(yvar), dataset)
  mn <- min(yall,na.rm=T)
  mx <- max(yall,na.rm=T)
  br <- seq(mn - 0.01*(mx-mn),mx + 0.01*(mx-mn),length=nbin)
  w <- br[2]-br[1]
  
  d <- split(dataset, dataset$pftlong)
  
  if(openwin)windows(4,7)
  par(mfrow=c(5,1), mar=c(0,0,0,0), oma=c(5,5,2,2))
  for(i in 1:5){
  
    x <- d[[i]]
    Y <- eval(substitute(yvar),x)
    Y <- Y[!is.na(Y)]
    
    h <- hist(Y, breaks=br, plot=FALSE)
  
    plot(br, br, ylim=c(0,max(h$counts)), axes=FALSE, type='n')
    for(j in 1:length(h$counts)){
      n <- h$counts[j]
      m <- h$mids[j]
      if(n == 0)next
      rect(xleft=m-w/2, xright=m+w/2, ybottom=0, ytop=n,  border=NA,col=col[i])
    }
    if(log)
        magaxis(side=1, unlog=1, tcl=-0.4)
    else
        axis(1)
    
    axis(2)
    
    legend("left", names(d)[i],fill=palette()[i], cex=0.7,bty='n')
    if(meanline)abline(v=mean(Y), lwd=2)
  }

mtext(side=2, line=3, text=ylab, outer=T)
mtext(side=1, line=3, text=xlab, outer=T)
}


palette(rainbow(5))
dat <- droplevels(subset(dataset2, pftlong %in% 
                           c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                             "EG-temperate")))

to.pdf({
histbypftlong(lalf_astbh, subset(dat, h.t>1.3), xlab=expression(A[L]/A[Sbh]~~(m^2~m^-2)))
histbypftlong(lmlf_astbh, subset(dat, h.t>1.3), xlab=expression(M[F]/A[Sbh]~~(m^2~m^-2)))
}, width=4, height=7, filename="output/figures/piperatio_bh_hist_bypftlong.pdf")


to.pdf({
histbypftlong(lalf_astba, dat, xlab=expression(A[L]/A[Sba]~~(m^2~m^-2)))
histbypftlong(lmlf_astba, dat, xlab=expression(M[F]/A[Sba]~~(m^2~m^-2)))
}, width=4, height=7, filename="output/figures/piperatio_basal_hist_bypftlong.pdf")


to.pdf(histbypftlong(lh.t, dat, meanline=F, ylab="Plant height (m)"),
       filename="output/figures/height_hist_bypftlong.pdf",
       width=4, height=7)


# Now aggregated by species-studyName
datag <- summaryBy(. ~ Group, data=dat, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                   id= ~ pftlong + pft)


to.pdf({
  histbypftlong(lalf_astbh, subset(datag, h.t>1.3), nbin=20, ylab="Nr of species", 
                xlab=expression(A[L]/A[Sbh]~~(m^2~m^-2)))
  histbypftlong(lmlf_astbh, subset(datag, h.t>1.3), nbin=20, ylab="Nr of species", 
                xlab=expression(M[F]/A[Sbh]~~(m^2~m^-2)))
}, width=4, height=7, filename="output/figures/piperatio_bh_hist_bypftlong_aggr.pdf")


to.pdf({
  histbypftlong(lalf_astba, datag, nbin=10, ylab="Nr of species", 
                xlab=expression(A[L]/A[Sba]~~(m^2~m^-2)))
  histbypftlong(lmlf_astba, datag, nbin=10, ylab="Nr of species", 
                xlab=expression(M[F]/A[Sba]~~(m^2~m^-2)))
}, width=4, height=7, filename="output/figures/piperatio_basal_hist_bypftlong_aggr.pdf")






# LMF scaling
lmfplot <- function(p, axes=T, cex=1, ...){
  
  dat <- dataset2[dataset2$pftlong == p,]
  
  with(dat, plot(log10(h.t), log10(m.st), axes=FALSE, pch=16,cex=cex,
                 xlab="Plant height (m)",
                 ylab="Woody or foliage mass (kg)",
                 col=make.transparent("brown"),
                 ...))
  with(dat, points(log10(h.t), log10(m.lf), pch=16,cex=cex,
                   col=make.transparent("forestgreen")))
  if(axes)magaxis(side=1:2, unlog=1:2)
  
}

to.pdf({
  
  for(p in unique(dataset2$pftlong)){
    lmfplot(p, main=p)
  }
  
}, filename="output/figures/mstmlf_pftlong_byht.pdf")



mstmlf_ht <- function(){
  #windows(4,9)
  par(mfrow=c(7,1), mar=c(0,0,0,0), oma=c(5,5,2,2))
  
  for(p in unique(dataset2$pftlong)){
    lmfplot(p, xlim=log10(c(0.05,105)), ylim=c(-6,6),axes=F, cex=0.5)
    magaxis(side=1, unlog=1, labels=FALSE)
    magaxis(side=2, unlog=2, labels=TRUE)
    box()
    legend("topleft", p, bty='n')
  }
  magaxis(side=1, unlog=1, labels=TRUE)
  mtext(side=1, text="Plant height (m)", line=3, outer=TRUE)
  mtext(side=2, text="Leaf or woody biomass (kg)", line=3, outer=TRUE)
}
to.pdf(mstmlf_ht(), width=4, height=9,
       filename="output/figures/mlfmst_byht_pftlong.pdf")






dat <- droplevels(subset(dataset2, h.t > 1.3 & pftlong %in% 
                           c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                             "EG-temperate")))

datag <- summaryBy(. ~ Group, data=dat, FUN=mean, na.rm=TRUE, keep.names=TRUE,
                   id= ~ pftlong + pft)

se <- function(x,...)sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))
datagg <- summaryBy(. ~ pft, data=datag, FUN=c(mean,se), na.rm=TRUE)


Cols <- c("blue","red","forestgreen")
palette(alpha(Cols,0.3))


windows(5,7)
par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
with(datag, plot(lsla, lmlf_astba2, pch=16, cex=0.8, col=pft,axes=FALSE))
magaxis(2, unlog=2)
magaxis(1, unlog=1, labels=F)
with(datagg, points(lsla.mean, lmlf_astba2.mean, col=Cols[pft], pch=15, cex=1.4))
with(datagg, arrows(x0=lsla.mean-2*lsla.se, x1=lsla.mean+2*lsla.se, 
                    y0=lmlf_astba2.mean,y1=lmlf_astba2.mean,code=3, length=0.04,angle=90,
                    col=Cols[pft]))
with(datagg, arrows(x0=lsla.mean, x1=lsla.mean, 
                         y0=lmlf_astba2.mean-2*lmlf_astba2.se,
                         y1=lmlf_astba2.mean+2*lmlf_astba2.se,code=3, length=0.04,angle=90,
                         col=Cols[pft]))
          
with(datag, plot(lsla, lalf_astba2, pch=16, cex=0.8, col=pft, axes=FALSE))
with(datagg, points(lsla.mean, lalf_astba2.mean, col=Cols[pft], pch=15, cex=1.4))
magaxis(1:2, unlog=1:2)
with(datagg, arrows(x0=lsla.mean-2*lsla.se, x1=lsla.mean+2*lsla.se, 
                    y0=lalf_astba2.mean,y1=lalf_astba2.mean,code=3, length=0.04,angle=90,
                    col=Cols[pft]))
with(datagg, arrows(x0=lsla.mean, x1=lsla.mean, 
                    y0=lalf_astba2.mean-2*lalf_astba2.se,
                    y1=lalf_astba2.mean+2*lalf_astba2.se,code=3, length=0.04,angle=90,
                    col=Cols[pft]))




# CV decrease (almost halved!)
cvm <- function(x)with(x, sd(10^lmlf_astbh, na.rm=T)/mean(10^lmlf_astbh, na.rm=T))
cva <- function(x)with(x, sd(10^lalf_astbh, na.rm=T)/mean(10^lalf_astbh, na.rm=T))

cvm(dataset2)
cva(dataset2)


sapply(split(dataset2, dataset2$pft), cvm)
sapply(split(dataset2, dataset2$pft), cva)

sapply(split(dataset2, dataset2$pftlong), cvm)
sapply(split(dataset2, dataset2$pftlong), cva)


sapply(split(datag, datag$pft), cvm)
sapply(split(datag, datag$pft), cva)


