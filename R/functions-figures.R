
alpha <- function (colour, alpha = NA) {
    col <- col2rgb(colour, TRUE)/255
    if (length(colour) != length(alpha)) {
        if (length(colour) > 1 && length(alpha) > 1) {
            stop("Only one of colour and alpha can be vectorised")
        }
        if (length(colour) > 1) {
            alpha <- rep(alpha, length.out = length(colour))
        }
        else if (length(alpha) > 1) {
            col <- col[, rep(1, length(alpha)), drop = FALSE]
        }
    }
    alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
    
    new_col <- rgb(col[1, ], col[2, ], col[3, ], alpha)
    new_col[is.na(colour)] <- NA
    new_col
}


# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.08, inset.x=inset, inset.y=inset, log.y=FALSE, log.x=FALSE,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])

  if(log.x)
    x <- 10^x

  if(log.y)
    y <- 10^y

  text(x,y,txt,...)
}

#' Means by some grouping variable g, accounting for random effect R.
mixmean <- function(yvar, g, data, R="Group"){

  data$P <- data[,g]
  data$Y <- data[,yvar]
  data$R <- data[,R]

  f <- lmer(Y ~ P - 1 + (1|R), data=data)

  ci <- suppressMessages(confint(f)[-c(1,2),])
  ci <- as.data.frame(cbind(10^fixef(f),10^ci))
  names(ci) <- c("y","lci","uci")
  ci[,g] <- gsub("P","", rownames(ci))
  rownames(ci) <- NULL

  df <- data.frame(unique(data$P))
  names(df) <- g
  df <- merge(df, ci, all=TRUE)

  lets <- as.data.frame(cld(glht(f, linfct=mcp(P="Tukey")))$mcletters$Letters)
  names(lets) <- "signifletters"
  lets[,g] <- rownames(lets)
  rownames(lets) <- NULL
  df <- merge(df, lets, all=TRUE)

  return(df)
}


#' Plot means by plant functional type
#' @description Optionally makes two panels. Adds CI based on \code{mixmean}.
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

  dat <- dataset

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

  # Multiple comparison letters.
  lets1 <- y1$signifletters
  if(!panel1only)lets2 <- y2$signifletters

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
  if(axis1){
    axis(1,labels=FALSE)
  }
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
    mtext(side=1, text=xlab, line=3, outer=T, las=0)
    mtext(side=2, at = 0.25, text=ylab2, line=3, outer=T, las=0)
    mtext(side=2, at = 0.75, text=ylab1, line=3, outer=T, las=0)
  }

  mtext(side=4, text=main, line=2, outer=TRUE)
  if(setpar)par(o)
}



#' Not a generic function! Only works for our case.
lsmeansGam <- function(obj, dataset, x){

  y <- list()
  for(i in 1:3){
    y[[i]] <- predict(obj[[i]],
                      data.frame(pft=levels(dataset$pft)[i],
                                 X=x),
                      se.fit=TRUE)
  }

  f <- function(z){
    p <- 10^c(z$fit, z$fit - 2*z$se.fit, z$fit + 2*z$se.fit)
    names(p) <- c("fit","lwr","upr")
    return(p)
  }
  tab <- sapply(y,f)
  colnames(tab) <- levels(dataset$pft)
  cis <- t(tab[2:3,])

  lets <- cld_generic(cis, rownames(cis))$Letters

  return(list(tab=as.data.frame(t(tab)),siglets=lets))

}



#' Plot prediction from gam
#' Used to take single prediction from multiple gams, plot as function of something.
plotGamPred <- function(obj,   # object returned by smoothplot, with pft as grouping (list of 3 gams)
                        dat,
                        xvar="llma",
                        xpred,
                        pftvar="pft",
                        siglets="bottom",
                        xlab="",
                        xlim=c(0,25),main="",
                        ylab=NULL,
                        ylim=NULL,
                        axis1=TRUE,
                        xaxislabels=TRUE,
                        panel.expr=NULL,
                        legend.text=NULL,
                        legend.where="topleft",
                        addlegend=FALSE,...){

  dat$P <- dat[,pftvar]

  # Calculate mixmeans
  if(is.character(xvar)){
    X <- mixmean(xvar,pftvar,dat)
  }

  z <- lsmeansGam(obj, dat, xpred)

  plot(X$y, z$tab$fit, xlim=xlim,axes=FALSE, pch=19, col=my_cols(), cex=1.3,
       ylim=ylim,xlab=xlab,ylab=ylab,
       panel.first={
         arrows(x0=X$lci, x1=X$uci, y0=z$tab$fit,
                y1=z$tab$fit,code=3,angle=90,length=0.025,col=my_cols())
         arrows(x0=X$y, x1=X$y, y0=z$tab$lwr,
                y1=z$tab$upr,code=3,angle=90,length=0.025,col=my_cols())
       })
  axis(2)
  if(axis1)axis(1,labels=xaxislabels)
  box()

  u <- par()$usr
  if(siglets == "bottom"){
    text(X$y, u[3] + 0.0*(u[4]-u[3]), z$siglets, pos=3, cex=0.9)
  }

  if(!is.null(panel.expr))eval(panel.expr)

  if(is.null(legend.text))legend.text <- rownames(X)
  if(addlegend){
    legend(legend.where, legend.text,pch=19,col=my_cols(),
           cex=legend.cex,pt.cex=1.2, bty='n')
  }

}



#' Plot least-square means
#' @description Plots least-square means, adds CI bars and significance letters.
#' @param y An object returned by \code{lmerTest::lsmeans},
#' @param x Either a numeric vector, or an object returned by \code{mixmean}
lsmeansPlot <- function(y,x, ylim=NULL, col=palette(), ...){


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
                  code=3,angle=90,length=0.025,col=col)
           arrows(x0=x$y, x1=x$y, y0=lci, y1=uci, angle=90, code=3,
                  length=0.025, col=col)
         }, pch=19, col=col,...)
    u <- par()$usr
    text(x$y, u[3] + 0.0*(u[4]-u[3]), lets, pos=3, cex=0.9)

  } else {

    plot(x, Y, ylim=ylim, col=col, pch=19,
         panel.first= arrows(x0=x, x1=x, y0=lci, y1=uci, angle=90, code=3,
                                       length=0.025, col=col),
         ...)
    u <- par()$usr
    text(x, u[3] + 0.0*(u[4]-u[3]), lets, pos=3, cex=0.9)


  }

}


#' Function for smoothplot. Probably not use otherwise.
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


#' Plot a generalized additive model
#' @param x Variable for X axis (unquoted)
#' @param y Variable for Y axis (unquoted)
#' @param data Dataframe containing x and y
#' @param kgam the \code{k} parameter for smooth terms in gam.
#' @param R An optional random effect (quoted)
#' @param log Whether to add log axes for x or y (but no transformations are done).
#' @param fitoneline Whether to fit only
smoothplot <- function(x,y,g=NULL,data,
                            fittype=c("gam","lm"),
                            kgam=4,
                            R=NULL,
                            randommethod=c("lmer","aggregate"),
                            axes=TRUE,
                            log="xy",
                            fitoneline=FALSE,
                            pointcols=NULL,
                            linecols=NULL,
                            xlab=NULL, ylab=NULL,
                            polycolor=alpha("lightgrey",0.7),
                            plotit=TRUE,
                            pch=16,
                            add=FALSE,
                            plotpoints=TRUE,
                            panel.first=NULL,
                            ...){

  if(add)axes <- FALSE
  fittype <- match.arg(fittype)
  randommethod <- match.arg(randommethod)

  if(!is.null(substitute(g))){
    data$G <- as.factor(eval(substitute(g),data))
  } else {
    fitoneline <- TRUE
    data$G <- 1
  }
  data$X <- eval(substitute(x),data)
  data$Y <- eval(substitute(y),data)
  data <- droplevels(data)

  data <- data[!is.na(data$X) & !is.na(data$Y) & !is.na(data$G),]

  if(is.null(pointcols))pointcols <- palette()
  if(is.null(linecols))linecols <- palette()

  if(is.null(xlab))xlab <- substitute(x)
  if(is.null(ylab))ylab <- substitute(y)

  # If randommethod = aggregate, average by group and fit simple gam.
  if(!is.null(R) && randommethod == "aggregate"){
    data$R <- data[,R]

    data <- summaryBy(. ~ R, FUN=mean, na.rm=TRUE, keep.names=TRUE, data=data,
                      id=~G)
    R <- NULL
  }


  if(!fitoneline){

    d <- split(data, data$G)

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

  if(plotit){
    
    if(!add){
      if(plotpoints){
        with(data, plot(X, Y, axes=FALSE, pch=pch, col=pointcols[G],
                      xlab=xlab, ylab=ylab, panel.first=panel.first, ...))
      } else {
        with(data, plot(X, Y, axes=FALSE, type='n', 
                        xlab=xlab, ylab=ylab, panel.first=panel.first,  ...))
      }
    } else {
      if(plotpoints){
        with(data, points(X, Y, pch=pch, col=pointcols[G],...))
      }
    }
    
    
    
    if(axes){
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
    }
    
    for(i in 1:length(fits)){

      if(fittype == "gam"){
        nd <- data.frame(X=seq(hran[[i]][1], hran[[i]][2], length=101))
        if(!inherits(fits[[i]], "try-error")){
          p <- predict(fits[[i]],nd,se.fit=TRUE)
          addpoly(nd$X, p$fit-2*p$se.fit, p$fit+2*p$se.fit, col=polycolor)
          lines(nd$X, p$fit, col=linecols[i], lwd=2)
        }
      }
      if(fittype == "lm"){
        pval <- summary(fits[[i]])$coefficients[2,4]
        LTY <- if(pval < 0.05)1 else 5
        predline(fits[[i]], col=linecols[i], lwd=2, lty=LTY, polycolor=polycolor)
      }
    }
  }
return(invisible(fits))
}


histbypft <- function(yvar, pftvar, dataset,
                      nbin=100,
                      plotwhat=c("counts","density"),
                      log=TRUE,
                      col=1:5,
                      meanlinecol="black",
                      cicol=alpha("grey",0.6),
                      xlab=NULL,
                      ylab="Nr. individuals",
                      ylim=NULL,
                      xlim=NULL,
                      legend.text=NULL,
                      legend.cex=1,
                      xaxis=NULL,
                      meanline=TRUE,
                      Means=NULL,
                      overlay=FALSE){

  if(is.null(Means))meanline <- FALSE
  yall <- eval(substitute(yvar), dataset)
  plotwhat <- match.arg(plotwhat)
  dataset$Group <- eval(substitute(pftvar), dataset)
  mn <- min(yall,na.rm=T)
  mx <- max(yall,na.rm=T)
  br <- seq(mn - 0.01*(mx-mn),mx + 0.01*(mx-mn),length=nbin)
  w <- br[2]-br[1]

  d <- split(dataset, dataset$Group)
  if(length(meanlinecol) == 1)meanlinecol <- rep(meanlinecol,length(d))
  if(is.null(legend.text))legend.text <- names(d)

  if(is.null(xaxis))xaxis <- 1:length(d)


  for(i in 1:length(d)){

    x <- d[[i]]
    Y <- eval(substitute(yvar),x)
    Y <- Y[!is.na(Y)]

    h <- hist(Y, breaks=br, plot=FALSE)
    if(is.null(ylim))
      ylim <- c(0,max(h[[plotwhat]]))
    else
      Ylim <- ylim
    
    if(!overlay || (overlay & i == 1)){
      plot(br, br, ylim=ylim, xlim=xlim, axes=FALSE, type='n')
    }
    for(j in 1:length(h[[plotwhat]])){
      n <- h[[plotwhat]][j]
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

    if(!overlay && meanline){
        rect(xleft=log10(Means$lci[i]), xright=log10(Means$uci[i]),
             ybottom=0, ytop=u[4], col=cicol, border=NA)
        segments(x0=log10(Means$y[i]), x1=log10(Means$y[i]),
                 y0=0, y1=u[4], col=meanlinecol[i])
    }
  }

  if(overlay && meanline){
    for(i in 1:length(d)){ 
      del <- 0.1 * -u[3]
      rect(xleft=log10(Means$lci[i]), xright=log10(Means$uci[i]),
           ybottom=u[3]+del, ytop=-del, col=cicol, border=NA)
      arrows(x0=log10(Means$y[i]), x1=log10(Means$y[i]),
               y0=u[3]+del, y1=-del, col=meanlinecol[i], length=0.05)
    }
  }
  
  mtext(side=2, line=3, text=ylab, outer=T, las=3)
  mtext(side=1, line=3, text=xlab, outer=T, las=1)
}


addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.8),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}
predline <- function(fit, from=NULL, to=NULL, 
                     polycolor=alpha("lightgrey",0.8),lty=1,...){

  if(is.null(from))from <- min(fit$model[,2], na.rm=TRUE)
  if(is.null(to))to <- max(fit$model[,2], na.rm=TRUE)

  newdat <- data.frame(X = seq(from,to, length=101))
  names(newdat)[1] <- names(coef(fit))[2]

  pred <- as.data.frame(predict(fit, newdat, se.fit=TRUE, interval="confidence")$fit)

  if(lty == -1){
    pval <- summary(fit)$coefficients[2,4]
    lty <- ifelse(pval < 0.05, 1, 5)
  }
  
  addpoly(newdat[[1]], pred$lwr, pred$upr, col=polycolor)
  ablinepiece(fit, from=from, to=to,lty=lty, ...)

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

#'@title Add a line to a plot
#'@description As \code{abline}, but with \code{from} and \code{to} arguments.
#'If a fitted linear regression model is used as asn argument, it uses the min and max values of the data used to fit the model.
#'@param a Intercept (optional)
#'@param b Slope (optional)
#'@param reg A fitted linear regression model (output of \code{\link{lm}}).
#'@param from Draw from this X value
#'@param to Draw to this x value
#'@param \dots Further parameters passed to \code{\link{segments}}
#'@export
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){

  # Borrowed from abline
  if (!is.null(reg)) a <- reg

  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)

    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }

  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)

}



# Simple axes function, bases on magaxis and code therein. 
log10axes <- function(side=1:2, logged=NULL, labels=TRUE){
  
  magaxis(side=side, unlog=side, labels=FALSE)
  
  loggedaxes <- c(FALSE,FALSE)
  if(!is.null(logged)){
    loggedaxes[logged] <- TRUE  
  }
  
  for(ii in side){
    lims <- sort(par("usr")[(1 + (ii-1)*2):(2 + (ii-1)*2)])
    sci.tick <- maglab(10^lims, n = 5, log = TRUE, exptext = TRUE, 
                      crunch = TRUE, logpretty = TRUE, usemultloc = FALSE, 
                      prettybase = 10, hersh = FALSE)
    
    if(labels){
      lab <- do.call(expression, lapply(log10(sci.tick$labat), function(i) bquote(10^.(i))))
      
      if(loggedaxes[ii]){
        axis(side = ii, at = sci.tick$labat, tick = FALSE, 
                labels = lab, mgp = par("mgp"))
      } else {
        axis(side = ii, at = log10(sci.tick$labat), tick = FALSE, 
             labels = lab, mgp = par("mgp"))
      } 
    }
  }
}

