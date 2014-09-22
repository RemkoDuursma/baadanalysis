
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