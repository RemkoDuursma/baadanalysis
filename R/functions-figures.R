
mixmean <- function(yvar, g, dataset=dat){
  
  dataset$P <- dataset[,g]
  dataset$Y <- dataset[,yvar]
  
  f <- lmer(Y ~ P - 1 + (1|Group), data=dataset)
  
  ci <- suppressMessages(confint(f)[-c(1,2),])
  ci <- as.data.frame(cbind(10^fixef(f),10^ci))
  rownames(ci) <- gsub("P","", rownames(ci))
  names(ci) <- c("y","lci","uci")
  return(ci)
}



meansbypft <- function(yvar1, yvar2, pftvar, 
                           xvar="lsla",
                           panel1.expr=NULL,
                           panel2.expr=NULL,
                           setpar=TRUE,
                           Cols=NULL,
                           legend.text=NULL,
                           legend.where="topright",
                           legend.cex=0.7,
                           xlab=expression(Specific~leaf~area~(m^2~kg^-1)), 
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
  
  sla <- mixmean(xvar,pftvar)
  mlfastbh <- mixmean(yvar1,pftvar)
  alfastbh <- mixmean(yvar2,pftvar)
  
  if(is.null(Cols))
    Cols <- rainbow(length(unique(dat$P)))
  
  # SHOULD be updated
  lmfit1 <- lm(mlfastbh$y ~ sla$y)
  lmfit2 <- lm(alfastbh$y ~ sla$y)
  
  # model fit for multiple comparisons
  f1 <- lmer(Y1 ~ P - 1 + (1|Group), data=dat)
  f2 <- lmer(Y2 ~ P - 1 + (1|Group), data=dat)
  
  # Multiple comparison letters.
  lets1 <- cld(glht(f1, linfct=mcp(P="Tukey")))$mcletters$Letters
  lets2 <- cld(glht(f2, linfct=mcp(P="Tukey")))$mcletters$Letters
  
  o <- par(no.readonly=TRUE)
  if(setpar)par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
  plot(sla$y, mlfastbh$y, xlim=xlim,axes=FALSE, pch=19, col=Cols, cex=1.3,
       ylim=ylim1,xlab=xlab,ylab=ylab1,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=mlfastbh$y, y1=mlfastbh$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=sla$y, x1=sla$y, y0=mlfastbh$lci, y1=mlfastbh$uci,code=3,angle=90,length=0.025,col=Cols)
         if(addtrend[1])ablinepiece(lmfit1)
       })
  axis(2)
  box()
  u <- par()$usr
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets1, pos=3, cex=0.9)
  axis(1,labels=FALSE)
  if(!is.null(panel1.expr))eval(panel1.expr)
  
  if(is.null(legend.text))legend.text <- rownames(sla)
  legend(legend.where, legend.text,pch=19,col=Cols, cex=legend.cex,pt.cex=1.2)
  plot(sla$y, alfastbh$y, xlim=xlim,pch=19, col=Cols, cex=1.3,
       ylim=ylim2,xlab=xlab, ylab=ylab2,
       panel.first={
         arrows(x0=sla$lci, x1=sla$uci, y0=alfastbh$y, y1=alfastbh$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=sla$y, x1=sla$y, y0=alfastbh$lci, y1=alfastbh$uci,code=3,angle=90,length=0.025,col=Cols)
         if(addtrend[2])ablinepiece(lmfit2)
       })
  u <- par()$usr
  if(!is.null(panel2.expr))eval(panel2.expr)
  text(sla$y, u[3] + 0.0*(u[4]-u[3]), lets2, pos=3, cex=0.9)
  mtext(side=1, text=xlab, line=3, outer=T)
  mtext(side=2, at = 0.25, text=ylab2, line=3, outer=T)
  mtext(side=2, at = 0.75, text=ylab1, line=3, outer=T)
  
  mtext(side=4, text=main, line=2, outer=TRUE)
  par(o)
}



smoothplotbypft <- function(x,y,data,pointcols=alpha(c("blue","red","forestgreen"),0.3),
                            linecols=c("deepskyblue3","red","chartreuse3"), 
                            xlab=NULL, ylab=NULL,
                            ...){
  
  data$pft <- as.factor(data$pft)
  data <- droplevels(data)
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  d <- split(data, data$pft)
  
  hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=4)))
  
  with(data, plot(X, Y, axes=FALSE, pch=16, col=pointcols[pft],
                  xlab=xlab, ylab=ylab, ...)) 
  magaxis(side=1:2, unlog=1:2)
  
  for(i in 1:length(d)){
    nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
    if(!inherits(fits[[i]], "try-error")){
      p <- predict(fits[[i]],nd,se.fit=TRUE)
      addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=alpha("lightgrey",0.7))
      lines(nd$X, p$fit, col=linecols[i], lwd=2)
    }
  }
}


histbypft <- function(yvar, pftvar, dataset, 
                      nbin=100, 
                      log=TRUE, 
                      col=1:5,
                      xlab=NULL, 
                      ylab="Nr. individuals", 
                      legend.text=NULL,
                      legend.cex=1,
                      xaxis=NULL,
                      meanline=TRUE){
  
  
  yall <- eval(substitute(yvar), dataset)
  dataset$Group <- eval(substitute(pftvar), dataset)
  mn <- min(yall,na.rm=T)
  mx <- max(yall,na.rm=T)
  br <- seq(mn - 0.01*(mx-mn),mx + 0.01*(mx-mn),length=nbin)
  w <- br[2]-br[1]
  
  d <- split(dataset, dataset$Group)
  if(is.null(legend.text))legend.text <- names(d)
  
  if(is.null(xaxis))xaxis <- 1:length(d)
  
  for(i in 1:length(d)){
    
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
    
    if(i %in% xaxis){
      if(log)
        magaxis(side=1, unlog=1, tcl=-0.4)
      else
        axis(1)
    }
    axis(2)
    
    legend("topleft", legend.text[i], cex=legend.cex,bty='n',text.font=2)
    if(meanline)abline(v=mean(Y), lwd=2)
  }
  
  mtext(side=2, line=3, text=ylab, outer=T, las=3)
  mtext(side=1, line=3, text=xlab, outer=T, las=1)
}

fitgam <- function(X,Y,dfr, k=-1, model=1){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  dfr <- droplevels(dfr)
  g<- gam(Y ~ s(X, k=k), data=dfr)
  
  return(g)
}
addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.8),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}


gamplotandpred <- function(dat, pftvar, yvar, 
                           xlab=NULL,
                           ylab=NULL, 
                           Hs = c(1,10), plotwhich=1:2, pointCols=NULL,
                           lineCols=NULL, vlines=TRUE,legend=TRUE){
  
  o <- par(no.readonly = TRUE)
  dat$PFT <- dat[,pftvar]
  dat$Y <- dat[,yvar]
  dat <- droplevels(subset(dat, !is.na(lh.t) & !is.na(Y)))
  
  d <- split(dat, dat$PFT)
  
  hran <- lapply(d, function(x)range(x$lh.t, na.rm=TRUE))
  fits <- lapply(d, function(x)fitgam("lh.t","Y",x, k=4))
  pred <- lapply(fits, function(x)predict(x, newdata=data.frame(X=log10(Hs)), 
                                          se.fit=TRUE))
  
  if(is.null(pointCols))
    pointCols <- rainbow(nlevels(dat$PFT))
  
  if(is.null(lineCols))
    lineCols <- pointCols
  
  if(is.null(xlab))xlab <- "lh.t"
  if(is.null(ylab))ylab <- yvar
  
  
  if(1 %in% plotwhich){
    with(dat, plot(lh.t, Y, axes=FALSE, pch=16, col=pointCols[PFT], ylab=ylab,xlab=xlab))
    magaxis(side=1:2, unlog=1:2)
    
    for(i in 1:length(d)){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      p <- predict(fits[[i]],nd,se.fit=TRUE)
      addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=alpha("lightgrey",0.7))
      lines(nd$X, p$fit, col=lineCols[i], lwd=2)
    }
    if(vlines)abline(v=log10(Hs), col="dimgrey", lty=5)
    if(legend)legend("bottomleft", levels(dat$PFT), fill=pointCols, cex=0.7)
    
    logY <- sapply(pred, "[[", "fit")
    se <- sapply(pred, "[[", "se.fit")
    Y_cil <- 10^(logY - 2*se)
    Y_ciu <- 10^(logY + 2*se)
    Y <- 10^logY
  }
  
  if(2 %in% plotwhich){
    par(mar=c(7,5,2,2), las=2)
    b <- barplot2(Y[1,], beside=T, plot.ci=TRUE, ci.l=Y_cil[1,], ci.u=Y_ciu[1,],
                  main=paste0("Height = ",Hs[1],"m"), col=pointCols, ylab=yvar)
    
    b <- barplot2(Y[2,], beside=T, plot.ci=TRUE, ci.l=Y_cil[2,], ci.u=Y_ciu[2,],
                  main=paste0("Height = ",Hs[2],"m"), col=pointCols, ylab=yvar)
  }
  
  df1 <- as.data.frame(Y)
  df1$H <- Hs
  df2 <- as.data.frame(se)
  df2$H <- Hs
  df1 <- melt(df1, id="H", value.name=yvar, variable.name=pftvar)
  df2 <- melt(df2, id="H", value.name=paste0(yvar,"_SE"), variable.name=pftvar)  
  df <- merge(df1, df2)
  par(o)
  return(invisible(df))
}

plotg <- function(g, which, ...){
  
  if(which == 1)
    which <- 1:3
  else
    which <- 4:6
  
  dataset2$llma <- log10(1/(10^dataset2$lsla))
  lma <- mixmean("llma", "pft", dataset2)
  g$y <- g[,3]
  g$lci <- g[,3] - 2*g[,4]
  g$uci <- g[,3] + 2*g[,4]
  
  plot(lma$y, g$y[which], col=Cols, cex=1.3, pch=19,
       axes=FALSE,
       panel.first={
         arrows(x0=lma$lci, x1=lma$uci, y0=g$y[which], y1=g$y[which],code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=lma$y, x1=lma$y, y0=g$lci[which], y1=g$uci[which],code=3,angle=90,length=0.025,col=Cols)
       }, ... )
  axis(1)
  axis(2)
  box()
  
}




