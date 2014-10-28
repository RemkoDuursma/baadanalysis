
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



meansbypft <- function(yvar1, yvar2=NULL, pftvar, 
                           xvar="llma",
                           dataset,
                           axis1=TRUE,
                           panel1.expr=NULL,
                           panel2.expr=NULL,
                           setpar=TRUE,
                           Cols=NULL,
                           panel1only=FALSE,
                          addlegend=TRUE,
                           legend.text=NULL,
                           legend.where="topright",
                           legend.cex=0.7,
                           siglets=c("bottom","symbol"),
                           xlab=expression(Specific~leaf~area~(m^2~kg^-1)), 
                           xlim=c(0,25),main="",
                           ylab1=NULL, ylab2=NULL, 
                           ylim1=NULL, ylim2=NULL){
    
  if(pftvar == "pftlong"){
    dat <- droplevels(subset(dataset, pftlong %in% 
                               c("DA-temperate","EA-temperate","EA-tropical","EG-boreal",
                                 "EG-temperate")))
  } else {
    dat <- dataset
  }
  
  siglets <- match.arg(siglets)
  if(siglets == "symbol")library(maptools)
  
  dat$Y1 <- dat[,yvar1]
  if(!panel1only)dat$Y2 <- dat[,yvar2]
  dat$P <- dat[,pftvar]  
  
  # If xvar is a character, calculate mixmeans, otherwise it is already a numeric vector.
  if(is.character(xvar)){
    X <- mixmean(xvar,pftvar,dat)
  } else {
    X <- list(y = xvar, lci=rep(NA,length(xvar)), uci=rep(NA,length(xvar)))
  }
  
  # mix means of Y variables.
  y1 <- mixmean(yvar1,pftvar,dat)
  if(!panel1only)y2 <- mixmean(yvar2,pftvar,dat)
  
  if(is.null(Cols))
    Cols <- rainbow(length(unique(dat$P)))
  
  # model fit for multiple comparisons
  f1 <- lmer(Y1 ~ P - 1 + (1|Group), data=dat)
  if(!panel1only)f2 <- lmer(Y2 ~ P - 1 + (1|Group), data=dat)
  
  # Multiple comparison letters.
  lets1 <- cld(glht(f1, linfct=mcp(P="Tukey")))$mcletters$Letters
  if(!panel1only)lets2 <- cld(glht(f2, linfct=mcp(P="Tukey")))$mcletters$Letters
  
  if(setpar){
    o <- par(no.readonly=TRUE)
    par(mfrow=c(2,1), oma=c(5,5,2,2), mar=c(0,0,0,0))
  }
  plot(X$y, y1$y, xlim=xlim,axes=FALSE, pch=19, col=Cols, cex=1.3,
       ylim=ylim1,xlab=xlab,ylab=ylab1,
       panel.first={
         arrows(x0=X$lci, x1=X$uci, y0=y1$y, 
                y1=y1$y,code=3,angle=90,length=0.025,col=Cols)
         arrows(x0=X$y, x1=X$y, y0=y1$lci, 
                y1=y1$uci,code=3,angle=90,length=0.025,col=Cols)
       })
  axis(2)
  if(axis1)axis(1,labels=FALSE)
  box()
  
  u <- par()$usr
  if(siglets == "bottom"){
    text(X$y, u[3] + 0.0*(u[4]-u[3]), lets1, pos=3, cex=0.9)
  }
  if(siglets == "symbol"){
    pointLabel(X$y, y1$y, lets1, cex=0.9)
  }
  
  if(!is.null(panel1.expr))eval(panel1.expr)
  
  if(is.null(legend.text))legend.text <- rownames(X)
  if(addlegend){
    legend(legend.where, legend.text,pch=19,col=Cols, 
         cex=legend.cex,pt.cex=1.2, bty='n')
  }
  
  if(!panel1only){
    plot(X$y, y2$y, xlim=xlim,pch=19, col=Cols, cex=1.3,
         ylim=ylim2,xlab=xlab, ylab=ylab2,axes=FALSE,
         panel.first={
           arrows(x0=X$lci, x1=X$uci, y0=y2$y, 
                  y1=y2$y,code=3,angle=90,length=0.025,col=Cols)
           arrows(x0=X$y, x1=X$y, y0=y2$lci, 
                  y1=y2$uci,code=3,angle=90,length=0.025,col=Cols)
         })
    axis(2)
    box()
    if(axis1)axis(1)
    u <- par()$usr
    if(!is.null(panel2.expr))eval(panel2.expr)
    if(siglets == "bottom"){
      text(X$y, u[3] + 0.0*(u[4]-u[3]), lets2, pos=3, cex=0.9)
    }
    if(siglets == "symbol"){
      pointLabel(X$y, y2$y, lets2, cex=0.9)
    }
    mtext(side=1, text=xlab, line=3, outer=T)
    mtext(side=2, at = 0.25, text=ylab2, line=3, outer=T)
    mtext(side=2, at = 0.75, text=ylab1, line=3, outer=T)
  }
  
  mtext(side=4, text=main, line=2, outer=TRUE)
  if(setpar)par(o)
}


lsmeansPlot <- function(y,x,lets,ylim=NULL,...){
  
  
  if(class(y)[1] != "lsmeans")
    stop("Need object returned by lmerTest::lsmeans")
  
  uci <- 10^y$lsmeans.table[["Upper CI"]]
  lci <- 10^y$lsmeans.table[["Lower CI"]]
  Y <- 10^y$lsmeans.table[["Estimate"]]
  
  if(is.null(ylim))ylim <- c(0, max(uci))
  
  numx <- is.list(x) & "lci" %in% names(x)
  
  lets <- cld.lsmeans(y)$Letters
  
  if(numx){
    plot(x$y, Y, 
         ylim=ylim, 
         panel.first={
           arrows(x0=x$lci, x1=x$uci, y0=Y, y1=Y,
                  code=3,angle=90,length=0.025,col=Cols)
           arrows(x0=x$y, x1=x$y, y0=lci, y1=uci, angle=90, code=3, 
                  length=0.025, col=Cols)
         }, pch=19, col=Cols,...)
    u <- par()$usr
    text(x$y, u[3] + 0.0*(u[4]-u[3]), lets, pos=3, cex=0.9)
    
  } else {
    
    plot(x, Y, ylim=ylim, col=Cols, pch=19, 
         panel.first= arrows(x0=x, x1=x, y0=lci, y1=uci, angle=90, code=3, 
                                       length=0.025, col=Cols),
         ...)
    u <- par()$usr
    text(x, u[3] + 0.0*(u[4]-u[3]), lets, pos=3, cex=0.9)

    
  }
  
  

}


smoothplotbypft <- function(x,y,data,pointcols=alpha(c("blue","red","forestgreen"),0.3),
                            fittype=c("gam","lm"),
                            kgam=4,
                            R=NULL,
                            log="xy",
                            fitoneline=FALSE,
                            linecols=c("deepskyblue3","red","chartreuse3"), 
                            xlab=NULL, ylab=NULL,
                            ...){
  
  fittype <- match.arg(fittype)
  data$pft <- as.factor(data$pft)
  data <- droplevels(data)
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  
  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$pft),]
  
  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)
  
  d <- split(data, data$pft)
  
  
  
  if(!fitoneline){
    if(fittype == "gam"){
      fits <- lapply(d, function(x)try(fitgam("X","Y",x, k=kgam, R=R)))
      if(!is.null(R))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- lapply(d, function(x)lm(Y ~ X, data=x))
    }
    hran <- lapply(d, function(x)range(x$X, na.rm=TRUE))
  } else {
    if(fittype == "gam"){
      fits <- list(fitgam("X","Y",data, k=kgam, R=R))
      if(!is.null(R))fits <- lapply(fits, "[[", "gam")
    } else {
      fits <- list(lm(Y ~ X, data=data))
    }
    hran <- list(range(data$X, na.rm=TRUE))
    
  }
  
  with(data, plot(X, Y, axes=FALSE, pch=16, col=pointcols[pft],
                  xlab=xlab, ylab=ylab, ...))
  
  if(log=="xy")magaxis(side=1:2, unlog=1:2)
  if(log=="x"){
    magaxis(side=1, unlog=1)
    axis(2)
    box()
  }
  if(log=="y"){
    magaxis(side=2, unlog=2)
    axis(1)
    box()
  }
  if(log==""){
    axis(1)
    axis(2)
    box()
  }
  
  for(i in 1:length(fits)){
    
    if(fittype == "gam"){
      nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
      if(!inherits(fits[[i]], "try-error")){
        p <- predict(fits[[i]],nd,se.fit=TRUE)
        addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=alpha("lightgrey",0.7))
        lines(nd$X, p$fit, col=linecols[i], lwd=2)
      }
    }
    if(fittype == "lm"){
      pval <- summary(fits[[i]])$coefficients[2,4]
      LTY <- if(pval < 0.05)1 else 5
      predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY)
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
                      meanline=TRUE,
                      Means=NULL){
  
  if(is.null(Means))meanline <- FALSE
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
    
    
    u <- par()$usr
    text(x=u[1], y=0.96*u[4], legend.text[i], cex=legend.cex,font=2,pos=4)
    
    if(meanline){
      
      rect(xleft=log10(Means$lci[i]), xright=log10(Means$uci[i]),
           ybottom=0, ytop=max(h$counts), col=alpha("grey",0.6), border=NA)
      segments(x0=log10(Means$y[i]), x1=log10(Means$y[i]), 
               y0=0, y1=max(h$counts))
      
    }
  }
  
  mtext(side=2, line=3, text=ylab, outer=T, las=3)
  mtext(side=1, line=3, text=xlab, outer=T, las=1)
}

fitgam <- function(X,Y,dfr, k=-1, R=NULL){
  dfr$Y <- dfr[,Y]
  dfr$X <- dfr[,X]
  if(!is.null(R)){
    dfr$R <- dfr[,R]
    model <- 2
  } else model <- 1
  dfr <- droplevels(dfr)
  
  
  if(model ==1){
    g <- gam(Y ~ s(X, k=k), data=dfr)
  }
  if(model ==2){
    g <- gamm(Y ~ s(X, k=k), random = list(R=~1), data=dfr)
  }
  
  return(g)
}
addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.8),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}
predline <- function(fit, from=NULL, to=NULL, ...){
  
  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)
  
  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]
  
  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)
  
  addpoly(newdat[[1]], pred$lwr, pred$upr)
  ablinepiece(fit, from=from, to=to, ...)
  
}

gamplotandpred <- function(dat, pftvar, yvar, 
                           setpar=TRUE,
                           xlab=NULL,
                           ylab=NULL, 
                           Hs = c(1,10), plotwhich=1:2, pointCols=NULL,
                           lineCols=NULL, vlines=TRUE,legend=TRUE){
  
  if(setpar)o <- par(no.readonly = TRUE)
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
  if(setpar)par(o)
  return(invisible(list(df=df, fits=fits)))
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



# Modified from plotrix::hexagon
Hexagon <- function (x, y, xdiam = 1, ydiam=xdiam, center=TRUE, col = NA, border = "black") {
  
  xx <- c(x, x, x + xdiam/2, x + xdiam, x + xdiam, 
          x + xdiam/2)
  yy <- c(y + ydiam * 0.125, y + ydiam * 
            0.875, y + ydiam * 1.125, y + ydiam * 0.875, y + 
            ydiam * 0.125, y - ydiam * 0.125)
  
  if(center){
    xx <- xx - xdiam/2
    yy <- yy - ydiam/2
  }
  
  polygon(xx, yy, col = col, border = border)
}

# Hacky method to get spacing between points.
getdiams <- function(cells){
  cel <- as.data.frame(cells)
  cel$xcount <- with(cel, ave(x, y, FUN=length))
  z <- subset(cel, xcount == max(xcount))
  xdiam <- median(diff(z$x))/2
  
  cel$ycount <- with(cel, ave(y, x, FUN=length))
  z <- subset(cel, ycount == max(ycount))
  ydiam <- median(diff(z$y))/2
  
  return(list(xdiam=xdiam, ydiam=ydiam))
}




